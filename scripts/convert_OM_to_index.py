#!/usr/bin/env python3
import argparse
import re
import sys
from urllib.parse import unquote

def parse_args():
    p = argparse.ArgumentParser(
        description="Map OperonMapper operon membership to locus tags using a GFF."
    )
    p.add_argument("-i", "--input", required=True, help="OperonMapper operons file")
    p.add_argument("-g", "--gff", required=True, help="Conjoined GFF file")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV (two columns: operon_number, locus_tag)")
    return p.parse_args()

def parse_gff_locus_maps(gff_path):
    """
    Build lookup maps (file order preserved per key):
      - gene -> [locus_tag...]
      - Name -> [locus_tag...]
      - gene_synonym -> [locus_tag...]
    Only entries carrying locus_tag are indexed.
    All attribute values are URL-decoded (e.g., %28 -> '(').
    """
    maps = {
        "gene": {},
        "Name": {},
        "gene_synonym": {}
    }
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]

            kv = {}
            for field in attrs.split(";"):
                if not field:
                    continue
                if "=" in field:
                    k, v = field.split("=", 1)
                    # URL-decode attribute values so '(' , ')' , quotes, etc. match raw IdGene text
                    kv[k] = unquote(v)

            locus = kv.get("locus_tag")
            if not locus:
                continue

            for key in ("gene", "Name", "gene_synonym"):
                val = kv.get(key)
                if val:
                    maps[key].setdefault(val, []).append(locus)

    return maps

PGAPTMP_RE = re.compile(r"(pgaptmp_\d+)", re.IGNORECASE)

# Allow ANY base chars, with optional trailing _N index (N>=1).
# IdGene examples handled: 'tet(M)', 'tet(M)_2', "aac(6')", "aph(3')-IIIa_2"
BASE_INDEX_ANY_RE = re.compile(r"^(?P<base>.+?)(?:_(?P<idx>\d+))?$")

def idgene_to_locus(idgene, maps):
    """
    Resolve IdGene -> locus_tag.
    1) If it contains 'pgaptmp_####', return that directly.
    2) Else interpret as BASE[_N] and pick the Nth (1-based) match.
       Priority: gene=BASE  -> Name=BASE -> gene_synonym=BASE (exact matches only).
    """
    m = PGAPTMP_RE.search(idgene)
    if m:
        return m.group(1)

    m2 = BASE_INDEX_ANY_RE.match(idgene)
    if not m2:
        return None
    base = m2.group("base")
    idx = m2.group("idx")
    n = int(idx) if idx else 1

    for key in ("gene", "Name", "gene_synonym"):
        candidates = maps[key].get(base)
        if candidates:
            if 1 <= n <= len(candidates):
                return candidates[n - 1]
            else:
                return None
    return None

def _first_nonempty_field(line):
    # Prefer tabs (OperonMapper output is TSV-like)
    fields = [f.strip() for f in line.rstrip("\n").split("\t")]
    fields = [f for f in fields if f != ""]
    if not fields:
        fields = [f for f in re.split(r"\s+", line.strip()) if f]
    return fields[0] if fields else ""

def parse_operon_file(operon_path):
    """
    Yields (operon_number, idgene_string) for each gene row.
    Handles:
      - leading tabs (empty first column)
      - header lines
      - NAME[_N] style IdGene
    """
    current_op = None
    with open(operon_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            if raw.lstrip().startswith("Operon") or raw.lstrip().startswith("IdGene"):
                continue

            stripped = raw.strip()
            if stripped.isdigit():
                current_op = int(stripped)
                continue

            if current_op is None:
                sys.stderr.write(f"[warn] Gene row encountered before any operon id: {raw}")
                continue

            idgene = _first_nonempty_field(raw)
            if not idgene or idgene in {"Operon", "IdGene", "Type"}:
                continue

            yield current_op, idgene

def main():
    args = parse_args()
    maps = parse_gff_locus_maps(args.gff)

    out_lines = []
    failures = 0
    seen_rows = 0

    for operon_num, idgene in parse_operon_file(args.input):
        seen_rows += 1
        locus = idgene_to_locus(idgene, maps)
        if locus is None:
            sys.stderr.write(f"[warn] Could not resolve locus_tag for IdGene '{idgene}' in operon {operon_num}\n")
            failures += 1
            continue
        out_lines.append(f"{operon_num}\t{locus}")

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("\n".join(out_lines) + ("\n" if out_lines else ""))

    sys.stderr.write(f"[info] Parsed {seen_rows} gene rows; wrote {len(out_lines)} lines; unresolved {failures}.\n")

if __name__ == "__main__":
    main()
