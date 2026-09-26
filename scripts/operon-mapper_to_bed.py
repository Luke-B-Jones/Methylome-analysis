#!/usr/bin/env python3
"""
operon2bed_with_gff_v2.py — Convert Operon Mapper output into BED-like operon intervals,
resolving contigs from a companion GFF3 with smarter ID handling and coordinate matching.

Output columns (5):
<contig>  <start>  <end>  <strand>  <operon_id>

Key features:
- Strips trailing _<num> from IdGene before GFF lookup (e.g., ltrA_2 -> ltrA).
- GFF is indexed on ID, Name, locus_tag, and gene; keeps ALL candidates.
- If multiple GFF hits match an identifier, picks by exact coordinate match, else by max overlap.
- BED semantics by default (0-based start, half-open end). Use --one-based to switch.

Mixed-contig handling: --on-mixed {dot,skip,first,majority} (default: dot)
"""

import argparse
import re
import sys
from collections import defaultdict, Counter
from urllib.parse import unquote

# -------------------- CLI --------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert Operon Mapper output into BED-like intervals with contigs from GFF3."
    )
    p.add_argument("-i", "--input", required=True, help="Operon Mapper output text file")
    p.add_argument("-o", "--output", required=True, help="Output BED-like file")
    p.add_argument("--gff", required=True, help="GFF3 file to resolve contigs (seqids)")
    p.add_argument("--one-based", action="store_true",
                   help="Write 1-based starts (DO NOT subtract 1). Default is 0-based start.")
    p.add_argument("--on-mixed", choices=["dot", "skip", "first", "majority"], default="dot",
                   help="If an operon spans multiple contigs: "
                        "'dot' -> write '.'; 'skip' -> omit; 'first' -> first seen; 'majority' -> majority contig.")
    p.add_argument("--require-contig", action="store_true",
                   help="Drop operons where any gene has unknown contig (instead of writing '.').")
    p.add_argument("--accept-types", default="CDS,gene,mRNA",
                   help="Comma-separated GFF types to index (default: CDS,gene,mRNA).")
    p.add_argument("--min-genes", type=int, default=1,
                   help="Minimum number of genes required in an operon to output it (default: 1).")
    return p.parse_args()

# -------------------- Operon-table parsing --------------------

def is_operon_header_line(line: str) -> bool:
    header_keywords = {"operon", "idgene", "type", "coggene", "posleft", "posright", "strand", "function"}
    toks = re.split(r"\s+", line.strip().lower())
    return sum(1 for t in toks if t in header_keywords) >= 3

def is_operon_id_line(line: str) -> bool:
    return re.fullmatch(r"\s*\d+\s*", line) is not None

def split_fields(line: str):
    s = line.strip()
    f = s.split("\t")
    if len(f) < 6:
        f = re.split(r"\s+", s)
    return f

def extract_gene_fields(fields):
    """
    Expect at least:
      0: IdGene
      1: Type
      2: COGgene
      3: PosLeft
      4: PosRight
      5: Strand
    Returns (idgene:str, left:int, right:int, strand:str) or None if malformed.
    """
    if len(fields) < 6:
        return None
    idgene = fields[0]
    try:
        left = int(fields[3])
        right = int(fields[4])
        strand = fields[5]
        if strand not in {"+", "-", "."}:
            ints = [int(x) for x in fields if re.fullmatch(r"[+-]?\d+", x)]
            if len(ints) >= 2:
                left, right = ints[0], ints[1]
            st = [x for x in fields if x in {"+", "-", "."}]
            strand = st[0] if st else "."
        if right < left:
            left, right = right, left
        return idgene, left, right, strand
    except ValueError:
        return None

# -------------------- GFF indexing and matching --------------------

def parse_gff_attributes(attr_str: str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = unquote(v)
    return attrs

def build_gff_index(gff_path: str, accept_types: set):
    """
    Build a multi-map: key -> list of candidate dicts with contig & coords.
    Keys include ID, Name, locus_tag, gene ( and any comma-split values ).
    Value item: {"contig":seqid, "start":int, "end":int, "strand":'+/-/.'}
    """
    idx = defaultdict(list)
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype not in accept_types:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
                if end_i < start_i:
                    start_i, end_i = end_i, start_i
            except ValueError:
                continue
            ad = parse_gff_attributes(attrs)
            keys = []
            for key in ("ID", "Name", "locus_tag", "gene"):
                if key in ad and ad[key]:
                    keys.extend(ad[key].split(","))
            rec = {"contig": seqid, "start": start_i, "end": end_i, "strand": strand}
            for k in keys:
                idx[k].append(rec)
    return idx

_UNDERSUFFIX = re.compile(r"^(.*?)(?:_\d+)$")

def strip_numeric_suffix(s: str) -> str:
    """
    ltrA_2 -> ltrA; leaves other IDs unchanged.
    """
    m = _UNDERSUFFIX.match(s)
    return m.group(1) if m else s

def interval_overlap(a_start, a_end, b_start, b_end):
    """
    Closed interval overlap length (GFF and operon table are 1-based inclusive).
    """
    lo = max(a_start, b_start)
    hi = min(a_end, b_end)
    return max(0, hi - lo + 1)

def choose_best_gff_hit(cands, left, right):
    """
    Pick the best candidate by:
    1) exact start/end match
    2) max overlap length
    3) first candidate (deterministic)
    """
    if not cands:
        return None
    # exact match first
    for rec in cands:
        if rec["start"] == left and rec["end"] == right:
            return rec
    # else max overlap
    best = None
    best_ov = -1
    for rec in cands:
        ov = interval_overlap(left, right, rec["start"], rec["end"])
        if ov > best_ov:
            best_ov = ov
            best = rec
    return best if best is not None else cands[0]

def resolve_gene_contig(idgene, left, right, gff_index):
    """
    Resolve contig for a gene entry using ID variations and coordinate disambiguation.
    Tries keys in this order for each name variant: exact, stripped prefix variants.
    """
    # 1) try exact IdGene
    cands = gff_index.get(idgene, [])
    if cands:
        hit = choose_best_gff_hit(cands, left, right)
        if hit:
            return hit["contig"]

    # 2) rescue common prefixes: cds- and gene-
    for prefix in ("cds-", "gene-"):
        k = idgene[len(prefix):] if idgene.startswith(prefix) else None
        if k:
            cands = gff_index.get(k, [])
            if cands:
                hit = choose_best_gff_hit(cands, left, right)
                if hit:
                    return hit["contig"]

    # 3) strip trailing _<num> (e.g., ltrA_2 -> ltrA)
    stripped = strip_numeric_suffix(idgene)
    if stripped != idgene:
        # Try stripped directly
        cands = gff_index.get(stripped, [])
        if cands:
            hit = choose_best_gff_hit(cands, left, right)
            if hit:
                return hit["contig"]
        # Try with common prefixes on stripped
        for prefix in ("cds-", "gene-"):
            k = stripped[len(prefix):] if stripped.startswith(prefix) else None
            if k:
                cands = gff_index.get(k, [])
                if cands:
                    hit = choose_best_gff_hit(cands, left, right)
                    if hit:
                        return hit["contig"]

    # 4) No match
    return None

# -------------------- Main --------------------

def main():
    args = parse_args()
    accept_types = set(x.strip() for x in args.accept_types.split(",") if x.strip())

    # Build GFF index: key -> [ {contig,start,end,strand}, ... ]
    gff_index = build_gff_index(args.gff, accept_types)

    # Parse operon table
    operon_genes = defaultdict(list)  # operon_id -> list of dicts
    current_operon = None
    saw_header = False

    with open(args.input, "r", encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            if not saw_header and is_operon_header_line(raw):
                saw_header = True
                continue
            if is_operon_id_line(raw):
                current_operon = int(raw.strip())
                continue
            if current_operon is None:
                continue
            fields = split_fields(raw)
            parsed = extract_gene_fields(fields)
            if parsed is None:
                continue
            idgene, left, right, strand = parsed
            operon_genes[current_operon].append(
                {"idgene": idgene, "left": left, "right": right, "strand": strand}
            )

    # Emit BED-like
    out = open(args.output, "w", encoding="utf-8")
    dropped = 0
    mixed_warned = 0
    unknown_warned = 0

    for operon_id in sorted(operon_genes.keys()):
        genes = operon_genes[operon_id]
        if len(genes) < args.min_genes:
            continue

        min_left = min(g["left"] for g in genes)
        max_right = max(g["right"] for g in genes)

        # Strand majority vote
        strands = [g["strand"] for g in genes if g["strand"] in {"+", "-"}]
        strand = "."
        if strands:
            c = Counter(strands)
            if c["+"] > c["-"]:
                strand = "+"
            elif c["-"] > c["+"]:
                strand = "-"

        # Resolve contigs per gene (with suffix stripping and coord disambiguation)
        contigs = []
        unknown_ids = []
        for g in genes:
            contig = resolve_gene_contig(g["idgene"], g["left"], g["right"], gff_index)
            if contig is None:
                unknown_ids.append(g["idgene"])
            else:
                contigs.append(contig)

        if unknown_ids:
            if args.require_contig:
                dropped += 1
                continue
            chosen_contig = "."
            if not unknown_warned:
                sys.stderr.write(
                    f"[warn] Some genes lacked contig mapping in GFF (showing '.' for contig). "
                    f"First examples: {', '.join(unknown_ids[:5])}\n"
                )
                unknown_warned = 1
        else:
            contig_counts = Counter(contigs)
            if len(contig_counts) == 1:
                chosen_contig = contigs[0]
            else:
                if args.on_mixed == "dot":
                    chosen_contig = "."
                elif args.on_mixed == "skip":
                    dropped += 1
                    continue
                elif args.on_mixed == "first":
                    chosen_contig = contigs[0]
                elif args.on_mixed == "majority":
                    chosen_contig = contig_counts.most_common(1)[0][0]
                else:
                    chosen_contig = "."
                if not mixed_warned:
                    sys.stderr.write(
                        "[warn] Some operons span multiple contigs; applying "
                        f"--on-mixed={args.on_mixed}. Example operon: {operon_id}\n"
                    )
                    mixed_warned = 1

        start = min_left if args.one_based else (min_left - 1)
        end = max_right
        out.write(f"{chosen_contig}\t{start}\t{end}\t{strand}\t{operon_id}\n")

    out.close()

    if dropped:
        sys.stderr.write(f"[info] Dropped {dropped} operon(s) by policy (mixed or unknown contig).\n")


if __name__ == "__main__":
    main()
