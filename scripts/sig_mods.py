#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
from collections import defaultdict

def open_maybe_gzip(path, mode="rt"):
    """Open plain text or .gz file transparently."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def iter_loci(path):
    """
    Yield (contig, start, end, strand) from a filtered per-base BED-like file.
    - Skips blank lines and comments.
    - Skips header/malformed lines by requiring numeric start/end.
    - Requires at least 6 columns to read strand; strand must be '+' or '-'.
    """
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            # Header / malformed guard: start/end must be integers
            try:
                start = int(parts[1])
                end   = int(parts[2])
            except ValueError:
                # header like: chrom  start  end ...
                continue
            # Strand (column 6 if present)
            strand = None
            if len(parts) >= 6:
                s = parts[5].strip()
                if s in {"+", "-"}:
                    strand = s
            if strand is None:
                # No valid strand column -> cannot participate in exact (contig,start,end,strand) match
                # Skip to avoid incorrect crossings.
                continue
            contig = parts[0].strip()
            yield (contig, start, end, strand)

def parse_args():
    p = argparse.ArgumentParser(
        description="Find per-base loci present in >1 treatment (exact match on contig,start,end,strand)."
    )
    p.add_argument(
        "-i", dest="inputs", nargs=2, action="append", metavar=("LABEL", "FILE"),
        required=True, help="Treatment label and input file (BED-like, header allowed)."
    )
    p.add_argument(
        "-o", dest="out", required=True,
        help="Output TSV (you may use .bed extension if preferred)."
    )
    p.add_argument(
        "--min-occ", type=int, default=2,
        help="Minimum number of distinct treatments a locus must appear in (default: 2)."
    )
    p.add_argument(
        "--no-header", action="store_true",
        help="Do not write a header row to the output."
    )
    return p.parse_args()

def main():
    args = parse_args()

    # Map locus -> set of labels where it appears
    locus2labels = defaultdict(set)

    # Also track per-treatment unique loci to avoid duplicate lines inflating counts
    per_label_sets = {}

    for label, path in args.inputs:
        if not os.path.exists(path):
            sys.stderr.write(f"[warn] Missing file for {label}: {path}\n")
            per_label_sets[label] = set()
            continue
        uniq = set()
        for locus in iter_loci(path):
            uniq.add(locus)
        per_label_sets[label] = uniq
        for locus in uniq:
            locus2labels[locus].add(label)

    # Filter loci by minimum distinct treatments
    qualifying = [(k, v) for k, v in locus2labels.items() if len(v) >= args.min_occ]

    # Sort deterministically: contig, start, end, strand, then treatments
    qualifying.sort(key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3], sorted(kv[1])))

    # Write output
    with open(args.out, "w") as out_f:
        if not args.no_header:
            out_f.write("contig\tstart\tend\tstrand\ttreatments\tcount\n")
        for (contig, start, end, strand), labels in qualifying:
            labs = ",".join(sorted(labels))
            out_f.write(f"{contig}\t{start}\t{end}\t{strand}\t{labs}\t{len(labels)}\n")

    # Minimal sanity report to stderr (does not pollute the TSV)
    sys.stderr.write("# Per-treatment unique locus counts:\n")
    for label, s in per_label_sets.items():
        sys.stderr.write(f"{label}\t{len(s)}\n")
    sys.stderr.write(f"[ok] Wrote {len(qualifying)} overlapping loci (≥{args.min_occ} treatments) to {args.out}\n")

if __name__ == "__main__":
    main()
