#!/usr/bin/env python3
"""
entropy_zero_pct.py

Recursively search a root directory for files named 'regions.bed'.
For each file:
- use column 5 (1-based; index 4) as mean entropy
- count rows where mean entropy is effectively 0.00
- compute % of rows (across all strands / all rows)

Output TSV:
sample_name    pct_effectively_zero

Definition of "effectively 0":
    abs(mean_entropy) < 0.005
which corresponds to displaying as 0.00 to 2 decimal places.
"""

import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Calculate per-file % of sites with mean entropy effectively 0.00."
    )
    p.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Root directory to recursively search for regions.bed files"
    )
    p.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )
    return p.parse_args()


def is_effectively_zero(entropy_value: float) -> bool:
    return abs(entropy_value) < 0.005


def analyze_regions_bed(path: Path):
    total_rows = 0
    zero_rows = 0

    with path.open("r", encoding="utf-8") as fh:
        for line_num, line in enumerate(fh, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                # Not enough columns to safely read mean entropy + strand
                continue

            try:
                mean_entropy = float(parts[4])   # 5th column (1-based)
            except ValueError:
                # Skip malformed numeric rows
                continue

            total_rows += 1
            if is_effectively_zero(mean_entropy):
                zero_rows += 1

    return total_rows, zero_rows


def main():
    args = parse_args()
    root = Path(args.input_dir)

    if not root.exists():
        raise FileNotFoundError(f"Input directory does not exist: {root}")

    regions_files = sorted(root.rglob("regions.bed"))

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("sample_name\tpct_effectively_zero\n")

        for bed_file in regions_files:
            sample_name = bed_file.parent.name
            total_rows, zero_rows = analyze_regions_bed(bed_file)

            if total_rows == 0:
                pct = "NA"
            else:
                pct = f"{(zero_rows / total_rows) * 100:.2f}"

            out.write(f"{sample_name}\t{pct}\n")


if __name__ == "__main__":
    main()