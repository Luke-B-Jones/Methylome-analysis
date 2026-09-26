#!/usr/bin/env python3
"""
operon_mapper_prep.py

Make BOTH outputs you asked for in one run:

1) **Operon-Mapper–ready GFF** (CDS-only, seqids matched to your genome FASTA;
   source column set to '.' unless you keep original).
2) **Per-locus DNA FASTA** extracted from the genome using the GFF
   (CDS-concatenated by default; respects strand and coordinates).

It also (optionally) writes a sanitized copy of the genome FASTA with short,
safe headers (<=35 chars, no spaces) so Operon-Mapper won't truncate/rename.

USAGE (matches your example; --protien is accepted as a no-op alias):
  python operon_mapper_prep.py \
    --protien /home/.../PGAP/annot.gff \
    --genome  "$GENOME" \
    --gff     /home/.../PGAP/annot.gff \
    --o-gff   /home/.../operon-mapper/operon-mapper_annot.gff \
    --o-nucl  /home/.../operon-mapper/DNA_annot.fasta

OPTIONAL (recommended): also write a sanitized genome FASTA and make GFF match it:
  --o-fasta /path/to/om_genome.fna

If your GFF uses placeholders (e.g., 'chromosome', 'plasmidA') that need to map
to real FASTA headers, provide explicit mappings (comma-separated old=new):
  --map 'chromosome=NZ_CP012345.1,plasmidA=NZ_CP012345.2'

Other options:
  --keep-non-cds           keep all GFF feature types (default: CDS only)
  --gff-source {.,orig}    set GFF source column; '.' (default) or 'orig'
  --feature-priority {CDS-first,gene-only}   how to extract per-locus DNA (default: CDS-first)
"""

import argparse
import sys
import re
from collections import defaultdict

# ===== FASTA utilities =====

SAFE_CHARS = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|")
MAX_HDR = 35
DNA_COMP = dict(
    A="T", C="G", G="C", T="A",
    R="Y", Y="R", S="S", W="W", K="M", M="K",
    B="V", D="H", H="D", V="B", N="N",
    a="t", c="g", g="c", t="a", r="y", y="r", s="s", w="w", k="m", m="k",
    b="v", d="h", h="d", v="b", n="n", **{"-":"-"}
)

def load_fasta_records(path):
    """Return list of (header, seq) where header is the first token after '>' and seq is uppercased."""
    recs = []
    h = None
    buf = []
    try:
        with open(path, "r", errors="replace") as fh:
            for ln in fh:
                if ln.startswith(">"):
                    if h is not None:
                        recs.append((h, "".join(buf).upper()))
                    h = ln[1:].split()[0]
                    buf = []
                else:
                    s = ln.strip()
                    if s:
                        buf.append(s)
            if h is not None:
                recs.append((h, "".join(buf).upper()))
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: FASTA not found: {path}\n"); sys.exit(2)
    if not recs:
        sys.stderr.write("ERROR: No sequences read from FASTA.\n"); sys.exit(2)
    return recs

def write_fasta(path, records, width=60):
    with open(path, "w") as out:
        for h, seq in records:
            out.write(f">{h}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def sanitize_header(h, idx, used):
    """Make a short, safe header for Operon-Mapper (<=35 chars, no spaces; unique)."""
    core = "".join(ch for ch in h if ch in SAFE_CHARS)
    if not core or len(core) > MAX_HDR:
        core = f"seq{idx+1}"
    base = core
    k = 1
    while core in used:
        k += 1
        suff = f"_{k}"
        core = base[:max(1, MAX_HDR - len(suff))] + suff
    used.add(core)
    return core

def revcomp(seq):
    return "".join(DNA_COMP.get(b, "N") for b in reversed(seq))

# ===== GFF utilities =====

def parse_attrs(attr_str):
    d = {}
    for chunk in attr_str.split(";"):
        chunk = chunk.strip()
        if chunk and "=" in chunk:
            k, v = chunk.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def read_gff_rows(path):
    rows = []
    try:
        with open(path, "r", errors="replace") as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                rows.append(parts)
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: GFF not found: {path}\n"); sys.exit(2)
    return rows

def parse_map(map_str):
    m = {}
    if not map_str:
        return m
    for pair in map_str.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            sys.stderr.write(f"Bad mapping '{pair}', expected old=new\n"); sys.exit(2)
        old, new = pair.split("=", 1)
        m[old] = new
    return m

# ===== Per-locus extraction =====

LOCUS_TOKEN = re.compile(r'([A-Za-z]+tmp_\d{3,})')  # e.g., pgaptmp_000123

def collect_locus_intervals(gff_rows, seqid_map_for_extraction):
    """
    Build per-locus structures for DNA extraction:
      locus_tag -> { 'seqid': str, 'strand': '+/-', 'gene_spans': [(s,e)], 'cds_parts': [(s,e),...] }
    Applies seqid_map_for_extraction to map GFF seqids to genome FASTA headers.
    """
    # Pass 1: gene/pseudogene -> geneID -> locus_tag
    geneid_to_locus = {}
    for p in gff_rows:
        seqid, source, ftype, start, end, score, strand, phase, attrs = p
        if ftype not in ("gene", "pseudogene"):
            continue
        a = parse_attrs(attrs)
        gid = a.get("ID")
        ltag = a.get("locus_tag")
        if gid and ltag:
            geneid_to_locus[gid] = ltag

    by_locus = defaultdict(lambda: {"seqid": None, "strand": None, "gene_spans": [], "cds_parts": []})

    for p in gff_rows:
        seqid, source, ftype, start, end, score, strand, phase, attrs = p
        a = parse_attrs(attrs)

        # Resolve locus_tag
        ltag = a.get("locus_tag")
        if not ltag and a.get("Parent"):
            for pid in a["Parent"].split(","):
                pid = pid.strip()
                if pid in geneid_to_locus:
                    ltag = geneid_to_locus[pid]
                    break
        if not ltag:
            for key in ("protein_id", "Name"):
                if key in a:
                    m = LOCUS_TOKEN.search(a[key])
                    if m:
                        ltag = m.group(1)
                        break
        if not ltag:
            continue

        # Map seqid for extraction (so it matches provided genome headers)
        mapped_seqid = seqid_map_for_extraction.get(seqid, seqid)

        rec = by_locus[ltag]
        if rec["seqid"] is None:
            rec["seqid"] = mapped_seqid
        if rec["strand"] is None and strand in ("+", "-"):
            rec["strand"] = strand

        try:
            s, e = int(start), int(end)
        except ValueError:
            continue
        if s > e:
            s, e = e, s

        if ftype == "CDS":
            rec["cds_parts"].append((s, e))
        elif ftype in ("gene", "pseudogene"):
            rec["gene_spans"].append((s, e))

    return by_locus

def extract_intervals(seq, intervals, strand):
    if not intervals:
        return ""
    intervals = sorted(intervals)
    L = len(seq)
    chunks = []
    for s, e in intervals:
        ss, ee = max(1, s), min(L, e)
        if ss <= ee:
            chunks.append(seq[ss - 1:ee])
    dna = "".join(chunks)
    return revcomp(dna) if strand == "-" else dna

# ===== Main pipeline =====

def main():
    ap = argparse.ArgumentParser(description="Produce OM-ready GFF and per-locus DNA FASTA in one go.")
    ap.add_argument("--protein", "--protien", dest="protein", help="[optional/ignored] path to protein FAA")
    ap.add_argument("--genome", required=True, help="Genome FASTA (FNA; contigs/chromosomes)")
    ap.add_argument("--gff", required=True, help="Annotation GFF")
    ap.add_argument("--o-gff", required=True, help="Output OM-ready GFF (CDS-only by default)")
    ap.add_argument("--o-nucl", required=True, help="Output per-locus DNA FASTA (extracted from genome)")
    ap.add_argument("--o-fasta", help="(Optional) write a sanitized genome FASTA whose headers are <=35 chars")
    ap.add_argument("--map", help="Comma list of old=new seqid mappings to apply to GFF BEFORE any sanitization")
    ap.add_argument("--keep-non-cds", action="store_true", help="Keep all GFF feature types (default: CDS only)")
    ap.add_argument("--gff-source", choices=[".", "orig"], default=".", help="Set GFF source column (default '.')")
    ap.add_argument("--feature-priority", choices=["CDS-first", "gene-only"], default="CDS-first",
                    help="Per-locus DNA extraction rule (default: CDS-first)")
    args = ap.parse_args()

    # Load genome FASTA
    genome_recs = load_fasta_records(args.genome)
    orig_headers = [h for h, _ in genome_recs]
    genome_dict_by_orig = {h: seq for h, seq in genome_recs}

    # Optional sanitized FASTA (recommended for OM)
    sanitized_hmap = {}  # original -> sanitized (if o-fasta given), else identity
    if args.o_fasta:
        used = set()
        for i, (h, _seq) in enumerate(genome_recs):
            sanitized_hmap[h] = sanitize_header(h, i, used)
        sanitized_records = [(sanitized_hmap[h], seq) for h, seq in genome_recs]
        write_fasta(args.o_fasta, sanitized_records)
        fa_names_for_gff = set(sanitized_hmap.values())
    else:
        sanitized_hmap = {h: h for h in orig_headers}
        fa_names_for_gff = set(orig_headers)

    # Read and optionally pre-map GFF seqids
    gff_rows = read_gff_rows(args.gff)
    user_map = parse_map(args.map)

    # ---- Part A: write OM-ready GFF (seqids matched to FASTA headers; CDS-only by default) ----
    kept_rows = []
    kept_seqids = set()
    dropped_rows = 0

    for parts in gff_rows:
        seqid, source, ftype = parts[0], parts[1], parts[2]

        # Apply explicit user mapping (e.g., chromosome -> NZ_CP...)
        seqid_mapped = user_map.get(seqid, seqid)

        # Then map to sanitized header (or identity when not sanitizing)
        if seqid_mapped in sanitized_hmap:
            new_seqid = sanitized_hmap[seqid_mapped]
        else:
            # If single-sequence genome, we can map everything to that one
            if len(fa_names_for_gff) == 1:
                new_seqid = next(iter(fa_names_for_gff))
            else:
                dropped_rows += 1
                continue

        if not args.keep_non_cds and ftype != "CDS":
            continue

        parts2 = parts[:]
        parts2[0] = new_seqid
        if args.gff_source == ".":
            parts2[1] = "."
        # else 'orig' -> keep original source
        kept_rows.append(parts2)
        kept_seqids.add(new_seqid)

    with open(args.o_gff, "w") as outg:
        outg.write("##gff-version 3\n")
        outg.write("##source=operon_mapper_prep\n")
        for p in kept_rows:
            outg.write("\t".join(p) + "\n")

    # ---- Part B: per-locus DNA FASTA (extracted from the ORIGINAL genome FASTA) ----
    # Build seqid mapping for extraction so GFF seqids can be looked up in genome_dict_by_orig
    # We use ONLY user_map here (not sanitized), because extraction uses the original genome headers.
    seqid_map_for_extraction = {k: v for k, v in user_map.items()}

    by_locus = collect_locus_intervals(gff_rows, seqid_map_for_extraction)

    out_count = 0
    skipped_no_seqid = 0
    skipped_missing_contig = 0
    skipped_empty = 0

    with open(args.o_nucl, "w") as outdna:
        for ltag, info in sorted(by_locus.items()):
            seqid = info["seqid"]
            strand = info["strand"] if info["strand"] in ("+", "-") else "+"
            if not seqid:
                skipped_no_seqid += 1
                continue
            if seqid not in genome_dict_by_orig:
                # If single-sequence genome, auto-map to it for extraction
                if len(genome_dict_by_orig) == 1:
                    seqid = next(iter(genome_dict_by_orig.keys()))
                else:
                    skipped_missing_contig += 1
                    continue
            if args.feature_priority == "CDS-first" and info["cds_parts"]:
                ivals = info["cds_parts"]
            elif info["gene_spans"]:
                gmin = min(s for s, _ in info["gene_spans"])
                gmax = max(e for _, e in info["gene_spans"])
                ivals = [(gmin, gmax)]
            else:
                skipped_empty += 1
                continue
            seq = extract_intervals(genome_dict_by_orig[seqid], ivals, strand)
            if not seq:
                skipped_empty += 1
                continue
            outdna.write(f">{ltag}\n")
            for i in range(0, len(seq), 60):
                outdna.write(seq[i:i+60] + "\n")
            out_count += 1

    # ---- Summary ----
    sys.stderr.write("=== SUMMARY ===\n")
    sys.stderr.write(f"Genome FASTA records: {len(genome_recs)}\n")
    if args.o_fasta:
        sys.stderr.write(f"Sanitized FASTA written: {args.o_fasta}\n")
    sys.stderr.write(f"GFF rows read: {len(gff_rows)}; written to OM-GFF: {len(kept_rows)}; dropped (seqid mismatch): {dropped_rows}\n")
    sys.stderr.write(f"Per-locus DNA FASTA written: {args.o_nucl} (records: {out_count})\n")
    if skipped_no_seqid:
        sys.stderr.write(f"  Skipped loci lacking seqid/strand in GFF: {skipped_no_seqid}\n")
    if skipped_missing_contig:
        sys.stderr.write(f"  Skipped loci with GFF seqid not found in genome FASTA: {skipped_missing_contig}\n")
    if skipped_empty:
        sys.stderr.write(f"  Skipped loci with no extractable intervals: {skipped_empty}\n")
    if kept_seqids.isdisjoint(fa_names_for_gff):
        sys.stderr.write("WARNING: No overlap between GFF seqids and FASTA headers after rewrite; check --map and headers.\n")

if __name__ == "__main__":
    main()
