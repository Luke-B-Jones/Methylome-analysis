#!/usr/bin/env python3
"""
Usage:
    python3 make_bed.py [-f] input_tblout > output.bed

Convert an Infernal cmsearch --tblout (default table) into BED format.

This script reads a cmsearch --tblout file where each non‐header line has at least 10 columns:
    0: target name (sequence/contig)
    1: target accession
    2: query name (model name)
    3: query accession
    4: mdl
    5: mdl from
    6: mdl to
    7: seq from
    8: seq to
    9: strand
    …plus additional columns

For each non‐header line it will:
  1. Extract chromosome/contig from column 0.
  2. Parse seq‐from (col 7), seq‐to (col 8) and strand (col 9).  If seq‐from > seq‐to, swap them so start ≤ end.
  3. Use query name (col 2) as the BED “name”.
  4. Deduplicate hits by (chrom, start, end, strand); for duplicates, concatenate query names with “;”.
  5. Print a BED line with six tab‐separated columns:
       chrom  start  end  name  score  strand
     (score is set to 0)

If the -f/--filter flag is given, any hit whose query name contains “tRNA” or “rRNA” is skipped.
"""

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert cmsearch --tblout to BED with optional rRNA/tRNA filtering."
    )
    parser.add_argument("-f", "--filter", action="store_true",
                        help="Skip lines where query name contains 'tRNA' or 'rRNA'.")
    parser.add_argument("input_tblout", help="cmsearch --tblout file")
    return parser.parse_args()

def main():
    args = parse_args()

    # Debug: report filter status
    print(f"[DEBUG] Filtering {'enabled' if args.filter else 'disabled'}.", file=sys.stderr)

    hits = {}  # key: (chrom, start, end, strand) -> list of query names

    with open(args.input_tblout, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split()
            if len(fields) < 10:
                continue  # not enough columns

            chrom = fields[0]
            query_name = fields[2]

            # apply optional filter
            if args.filter and ('tRNA' in query_name or 'rRNA' in query_name):
                continue

            # parse coordinates
            try:
                s = int(fields[7])
                e = int(fields[8])
                strand = fields[9]
            except ValueError:
                continue  # malformed coords

            # ensure start <= end
            start, end = (s, e) if s <= e else (e, s)

            key = (chrom, start, end, strand)
            hits.setdefault(key, []).append(query_name)

    # output in sorted order
    for chrom, start, end, strand in sorted(hits.keys(), key=lambda k: (k[0], k[1], k[2], k[3])):
        names = ";".join(hits[(chrom, start, end, strand)])
        # BED: chrom, start, end, name, score, strand
        print(f"{chrom}\t{start}\t{end}\t{names}\t0\t{strand}")

if __name__ == "__main__":
    main()
