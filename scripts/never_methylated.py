#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

TARGET_MODS = {"4mC", "5mC", "6mA"}
TARGET_FILENAME = "control_treatment.bed"
ZERO_THRESHOLD = 0.005  # abs(value) < 0.005 counts as effectively 0.00


def find_output_label(file_path: Path, input_dir: Path) -> str | None:
    rel = file_path.relative_to(input_dir)
    parent_parts = rel.parts[:-1]  # exclude filename

    for i, part in enumerate(parent_parts):
        if part in TARGET_MODS:
            return Path(*parent_parts[i:]).as_posix()

    return None


def parse_number(value: str) -> float:
    """
    Parse either:
      5
      0.0049
      21839:2
      21839:0.67

    by taking the part after the last ':' when present.
    """
    s = str(value).strip()
    if not s:
        raise ValueError("empty value")

    if ":" in s:
        s = s.rsplit(":", 1)[1]

    return float(s)


def analyse_bed_file(file_path: Path) -> tuple[int, int]:
    """
    Returns:
      eligible_rows, zero_like_rows

    eligible_rows: rows where treatment_counts >= 5
    zero_like_rows: eligible rows where abs(control_pct_modified) < 0.005
    """
    eligible_rows = 0
    zero_like_rows = 0

    with file_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid file: {file_path}")

        reader.fieldnames = [h.strip() if h is not None else h for h in reader.fieldnames]

        required = {"treatment_counts", "control_pct_modified"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns in {file_path}: {', '.join(sorted(missing))}"
            )

        for line_no, row in enumerate(reader, start=2):
            try:
                treatment_counts = parse_number(row["treatment_counts"])
                control_pct_modified = parse_number(row["control_pct_modified"])
            except (TypeError, ValueError):
                print(
                    f"Warning: skipping non-numeric row in {file_path} at line {line_no}",
                    file=sys.stderr,
                )
                continue

            if treatment_counts >= 5:
                eligible_rows += 1
                if abs(control_pct_modified) < ZERO_THRESHOLD:
                    zero_like_rows += 1

    return eligible_rows, zero_like_rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively find control_treatment.bed files under 4mC/5mC/6mA folders, "
            "then calculate the percentage of rows with treatment_counts >= 5 where "
            "control_pct_modified is effectively 0.00."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input root directory")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    input_dir = Path(args.input).expanduser().resolve()
    output_file = Path(args.output).expanduser()

    if not input_dir.is_dir():
        parser.error(f"Input path is not a directory: {input_dir}")

    results: list[tuple[str, str]] = []

    for bed_file in sorted(input_dir.rglob(TARGET_FILENAME)):
        label = find_output_label(bed_file, input_dir)
        if label is None:
            continue

        eligible_rows, zero_like_rows = analyse_bed_file(bed_file)

        percentage_str = "NA"
        if eligible_rows > 0:
            percentage = (zero_like_rows / eligible_rows) * 100
            percentage_str = f"{percentage:.3f}%"

        results.append((label, percentage_str))

    with output_file.open("w", encoding="utf-8", newline="") as out:
        for label, percentage_str in sorted(results):
            out.write(f"{label}\t{percentage_str}\n")


if __name__ == "__main__":
    main()