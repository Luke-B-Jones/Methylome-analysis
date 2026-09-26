#!/usr/bin/env python3
import os, sys, argparse, gzip
from typing import Dict, List, Tuple, Iterator

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)

def detect_delim(sample: str) -> str:
    if "\t" in sample: return "\t"
    if "," in sample:  return ","
    return "\t"

def read_locus_list(table_path: str, col_index: int = 1) -> List[str]:
    """
    Read non-header rows from table, return values from 1-based column `col_index`.
    Skips header if first field is 'Locus_tag' (case-insensitive).
    """
    loci: List[str] = []
    with open(table_path, "r", encoding="utf-8-sig") as fh:
        first = fh.readline()
        if not first:
            return loci
        delim = detect_delim(first)
        fh.seek(0)
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split(delim)
            if i == 0 and parts and parts[0].strip().lower() == "locus_tag":
                continue  # header
            if len(parts) < col_index:
                continue
            tag = parts[col_index - 1].strip()
            if tag:
                loci.append(tag)

    # Deduplicate while preserving order
    seen = set()
    out: List[str] = []
    for x in loci:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def fasta_iter(handle) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_gt, sequence_string) for each FASTA record."""
    hdr = None
    seq_chunks: List[str] = []
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq_chunks)
            hdr, seq_chunks = line[1:].strip(), []
        else:
            seq_chunks.append(line)
    if hdr is not None:
        yield hdr, "".join(seq_chunks)

def index_fasta_tokens(fasta_path: str) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Map common header tokens -> sequence.
    - Primary key = first whitespace-delimited token (strip 'lcl|' prefix).
    - Also index pieces split by '|' for convenience.
    Returns (token_map, full_list) for optional substring fallback.
    """
    token_map: Dict[str, str] = {}
    full_list: List[Tuple[str, str]] = []
    with open_maybe_gzip(fasta_path, "rt") as fh:
        for hdr, seq in fasta_iter(fh):
            full_list.append((hdr, seq))
            token = hdr.split()[0] if hdr else ""
            if token.startswith("lcl|"):
                token = token[4:]
            if token and token not in token_map:
                token_map[token] = seq
            if "|" in token:
                for piece in token.split("|"):
                    p = piece.strip()
                    if p and p not in token_map:
                        token_map[p] = seq
    return token_map, full_list

def wrap(seq: str, width: int = 60) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def open_out(path: str):
    return gzip.open(path, "wt") if path.lower().endswith(".gz") else open(path, "w")

def main():
    ap = argparse.ArgumentParser(
        description="Create a single multi-FASTA for BLASTN: take locus IDs from -i (col 1) and extract matching sequences from -g."
    )
    ap.add_argument("-i", "--input-table", required=True,
                    help="TSV/CSV with Locus_tag in column 1 (header 'Locus_tag' optional).")
    ap.add_argument("-g", "--genome-fasta", required=True,
                    help="Multi-FASTA (can be .gz) containing sequences; headers’ first token should match locus IDs.")
    ap.add_argument("-o", "--output", required=True,
                    help="Output multi-FASTA file (can be .gz).")
    ap.add_argument("--col-index", type=int, default=1,
                    help="1-based column index for locus tag in the table (default: 1).")
    ap.add_argument("--wrap", type=int, default=60,
                    help="FASTA line wrap width (default: 60; set 0 for no wrap).")
    ap.add_argument("--no-substring-fallback", action="store_true",
                    help="Only use exact token matches; disable unique substring fallback.")
    args = ap.parse_args()

    loci = read_locus_list(args.input_table, args.col_index)
    if not loci:
        print("[ERROR] No locus tags found in input table.", file=sys.stderr)
        sys.exit(1)

    token_map, full_list = index_fasta_tokens(args.genome_fasta)
    headers = [h for h, _ in full_list]

    n_written = 0
    missing: List[str] = []

    with open_out(args.output) as out_f:
        for locus in loci:
            seq = token_map.get(locus)
            if seq is None and not args.no_substring_fallback:
                # Fallback: unique header containing the locus string
                idxs = [i for i, h in enumerate(headers) if locus in h]
                if len(idxs) == 1:
                    seq = full_list[idxs[0]][1]
            if seq is None:
                missing.append(locus)
                continue

            out_f.write(f">{locus}\n{wrap(seq, args.wrap)}\n")
            n_written += 1

    # Brief summary to stderr
    print(f"[DONE] Wrote {n_written} records to {args.output}. Missing: {len(missing)}", file=sys.stderr)
    if missing:
        print("[WARN] Missing IDs (first 20): " + ", ".join(missing[:20]), file=sys.stderr)

if __name__ == "__main__":
    main()
