#!/usr/bin/env python3
import os, sys, argparse, gzip, re
from typing import Iterator, Tuple, List, Dict

# Accept .fa/.fas/.fna/.fasta with optional .gz (case-insensitive)
FA_PAT = re.compile(r'\.(?:fa|fas|fna|fasta)(?:\.gz)?$', re.IGNORECASE)

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)

def detect_delim(line: str) -> str:
    if "\t" in line: return "\t"
    if "," in line:  return ","
    return "\t"

def fasta_iter(handle) -> Iterator[Tuple[str, str]]:
    hdr, seq = None, []
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line[1:].strip(), []
        else:
            seq.append(line)
    if hdr is not None:
        yield hdr, "".join(seq)

def parse_meta_ids_with_year(meta_path: str, id_col_idx: int = 1, year_col_idx: int = 3) -> List[str]:
    """
    Return IDs (strings) where year column is non-empty.
    Columns are 1-based indices. Header is optional and skipped if it looks like one.
    """
    ids: List[str] = []
    with open(meta_path, "r", encoding="utf-8-sig") as fh:
        first = fh.readline()
        if not first:
            return ids
        delim = detect_delim(first)
        fh.seek(0)
        lines = fh.readlines()

    start = 1 if ("id" in lines[0].lower() and "year" in lines[0].lower()) else 0
    i_idx, y_idx = id_col_idx - 1, year_col_idx - 1

    for line in lines[start:]:
        parts = line.rstrip("\n").split(delim)
        if len(parts) <= max(i_idx, y_idx):
            continue
        rid = parts[i_idx].strip()
        year = parts[y_idx].strip()
        if rid and year != "":
            ids.append(rid)
    return ids

def find_fasta_for_id(root: str, rid: str) -> List[str]:
    """Recursively find files whose basename starts with the ID and has a FASTA extension."""
    hits: List[str] = []
    for dp, _ds, files in os.walk(root):
        for fn in files:
            if fn.startswith(rid) and FA_PAT.search(fn):
                hits.append(os.path.join(dp, fn))
    return sorted(hits)

def count_contigs_only(fasta_path: str) -> int:
    """Fast pass: count lines that start with '>'."""
    n = 0
    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n

def write_contigs_renamed(fasta_path: str, genome_id: str, out_handle, start_idx: int = 0) -> Tuple[int, int]:
    """
    Write contigs from fasta_path as >{genome_id}_{k} where k continues from start_idx+1.
    Returns (num_contigs_written, total_bp_written).
    """
    n = 0
    bp = 0
    with open_maybe_gzip(fasta_path, "rt") as fh:
        for _hdr, seq in fasta_iter(fh):
            n += 1
            bp += len(seq)
            out_handle.write(f">{genome_id}_{start_idx + n}\n{seq}\n")
    return n, bp

def main():
    ap = argparse.ArgumentParser(
        description="Build a mega FASTA: keep IDs from metadata with non-empty year (col 3). "
                    "For each ID, find FASTA(s) under -i whose basename starts with the ID; "
                    "skip FASTAs with > --max contigs if provided; write headers as ID_NUM."
    )
    ap.add_argument("-i", "--input-root", required=True,
                    help="Root folder containing per-assembly FASTAs (searched recursively).")
    ap.add_argument("-t", "--metadata", required=True,
                    help="Metadata TSV/CSV. Column 1 = ID, Column 3 = year. Header optional.")
    ap.add_argument("-o", "--output", required=True,
                    help="Output FASTA file (e.g., /path/to/database.fasta).")
    ap.add_argument("--id-col-index", type=int, default=1,
                    help="1-based column index for ID in metadata (default: 1).")
    ap.add_argument("--year-col-index", type=int, default=3,
                    help="1-based column index for year in metadata (default: 3).")
    ap.add_argument("--first-only", action="store_true",
                    help="Use only the first matching FASTA per ID (skip additional matches).")
    ap.add_argument("--max", dest="max_contigs", type=int, default=None,
                    help="Skip any FASTA with more than this number of contigs (e.g., --max 100).")
    args = ap.parse_args()

    # Ensure output dir exists
    out_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    os.makedirs(out_dir, exist_ok=True)

    ids = parse_meta_ids_with_year(args.metadata, args.id_col_index, args.year_col_index)
    # Deduplicate while preserving order
    seen = set()
    ordered_ids = [x for x in ids if not (x in seen or seen.add(x))]

    total_ids_used = 0

    with open(args.output, "w") as out_f:
        for rid in ordered_ids:
            hits = find_fasta_for_id(args.input_root, rid)
            if not hits:
                continue

            if args.first_only and len(hits) > 1:
                hits = hits[:1]

            wrote_any_for_id = False
            next_index_for_id = 0

            for fp in hits:
                # Apply --max contigs filter if requested
                if args.max_contigs is not None:
                    c = count_contigs_only(fp)
                    if c > args.max_contigs:
                        continue

                n_written, _bp = write_contigs_renamed(fp, rid, out_f, start_idx=next_index_for_id)
                if n_written > 0:
                    wrote_any_for_id = True
                    next_index_for_id += n_written

            if wrote_any_for_id:
                total_ids_used += 1

    # ONLY print genomes count (stdout), as requested
    print(f"{total_ids_used}.")

if __name__ == "__main__":
    main()
