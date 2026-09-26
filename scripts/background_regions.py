#!/usr/bin/env python3
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="""Find regions in the genome not covered by any of the provided intervals.

Inputs (all BED files are 0-based half-open; e.g., "chrom 0 10 ..." covers bases 1..10 in 1-based inclusive):
  -p/--operon      : Operon BED file (0-based half-open; robust to extra columns; strand auto-detected).
  -g/--gene        : Gene BED file (0-based half-open; robust to extra columns; strand auto-detected).
  -u/--upstream    : Upstream regions BED file (0-based half-open; robust to extra columns; strand auto-detected).
  -3/--three-prime : (Optional) Terminator/3' regions BED file (0-based half-open; robust to extra columns; strand auto-detected).
  -f/--fai-file    : FAI file specifying contig sizes.
  -o/--out-bed     : Output BED file (0-based half-open) of uncovered regions.

Notes:
- This script merges coverage per-contig and per-strand (strand taken from any '+'/'-' column if present; defaults to '+').
- All inputs are treated as 0-based half-open and converted to 1-based inclusive for merging.
- Output intervals are in standard 0-based half-open BED with a strand column.
"""
    )
    parser.add_argument("-p", "--operon", required=True,
                        help="Operon BED file (0-based half-open).")
    parser.add_argument("-g", "--gene", required=True,
                        help="Genes BED file (0-based half-open).")
    parser.add_argument("-u", "--upstream", required=True,
                        help="Upstream regions BED file (0-based half-open).")
    parser.add_argument("-3", "--three_prime", required=False,
                        help="(Optional) Terminator/3' BED file (0-based half-open).")
    parser.add_argument("-f", "--fai_file", required=True,
                        help="FAI file specifying contig sizes.")
    parser.add_argument("-o", "--out_bed", required=True,
                        help="Output BED file (0-based half-open) of uncovered regions.")
    return parser.parse_args()

def detect_strand(parts):
    """
    Return '+' or '-' if present in any column (from col 3 onward), else '+'.
    """
    for tok in parts[3:]:
        if tok == '+' or tok == '-':
            return tok
    return '+'

def parse_bed0_to_1inclusive(bed_file):
    """
    Parse a BED-like file where coordinates are 0-based half-open.
    Convert to 1-based inclusive intervals and group by chrom and strand.
    Accepts extra columns; strand is auto-detected from any column as '+' or '-'.
    """
    intervals = defaultdict(lambda: defaultdict(list))
    with open(bed_file, 'r') as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                chrom = parts[0]
                start0 = int(parts[1])
                end0   = int(parts[2])
            except ValueError:
                # Malformed numeric fields; skip
                continue
            if end0 <= start0:
                # zero-length or invalid; skip
                continue
            # Convert 0-based half-open -> 1-based inclusive
            start1 = start0 + 1
            end1   = end0
            strand = detect_strand(parts)
            intervals[chrom][strand].append((start1, end1))
    return intervals

def parse_fai(fai_file):
    """
    Parse FAI file to map contig -> size.
    """
    sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom  = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                continue
            sizes[chrom] = length
    return sizes

def merge_intervals_1based(intervals):
    """
    Merge overlapping/adjacent intervals in 1-based inclusive space.
    """
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(intervals[0])]
    for current in intervals[1:]:
        last = merged[-1]
        # If current starts at or before last.end + 1, they touch/overlap
        if current[0] <= last[1] + 1:
            if current[1] > last[1]:
                last[1] = current[1]
        else:
            merged.append(list(current))
    return merged

def subtract_intervals_1based(contig_size, merged_intervals):
    """
    Given contig size and merged covered intervals (1-based inclusive),
    return gaps (uncovered regions) as 1-based inclusive intervals.
    """
    uncovered = []
    current = 1
    for (start, end) in merged_intervals:
        if start > current:
            uncovered.append((current, start - 1))
        current = end + 1
    if current <= contig_size:
        uncovered.append((current, contig_size))
    return uncovered

def main():
    args = parse_args()

    # Parse required sources (0-based half-open -> 1-based inclusive)
    operon   = parse_bed0_to_1inclusive(args.operon)
    genes    = parse_bed0_to_1inclusive(args.gene)
    upstream = parse_bed0_to_1inclusive(args.upstream)

    # Optional three_prime (now also 0-based half-open)
    three_prime = defaultdict(lambda: defaultdict(list))
    if args.three_prime:
        three_prime = parse_bed0_to_1inclusive(args.three_prime)

    contig_sizes = parse_fai(args.fai_file)

    # Accumulate coverage
    coverage = defaultdict(lambda: defaultdict(list))
    def add_coverage(source):
        for chrom, strand_dict in source.items():
            for strand, interval_list in strand_dict.items():
                coverage[chrom][strand].extend(interval_list)

    add_coverage(operon)
    add_coverage(genes)
    add_coverage(upstream)
    add_coverage(three_prime)

    # Compute uncovered regions per contig and per strand
    uncovered = defaultdict(lambda: defaultdict(list))
    for chrom, size in contig_sizes.items():
        for strand in ['+', '-']:
            cov_list = coverage[chrom][strand]
            if not cov_list:
                uncovered[chrom][strand].append((1, size))
                continue
            merged = merge_intervals_1based(cov_list)
            gaps = subtract_intervals_1based(size, merged)
            uncovered[chrom][strand] = gaps

    # Write output as 0-based half-open BED with strand
    # Convert 1-based inclusive -> 0-based half-open: start0 = start1 - 1; end0 = end1
    with open(args.out_bed, 'w') as out:
        for chrom, strand_dict in uncovered.items():
            for strand, intervals_list in strand_dict.items():
                for i, (start1, end1) in enumerate(intervals_list, 1):
                    bed_start = start1 - 1
                    bed_end   = end1
                    name = f"uncovered_{chrom}_{strand}_{i}"
                    out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}\n")

    print(f"[INFO] Wrote uncovered regions to {args.out_bed}")

if __name__ == "__main__":
    main()
