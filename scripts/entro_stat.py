#!/usr/bin/env python3
from __future__ import annotations
import sys, argparse, gzip
from typing import Dict, List, Tuple, Optional

def smart_open(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def _is_strand_token(tok: str) -> bool:
    return tok in {"+", "-"}

def load_bed_name_to_strand(
    bed_path: str,
    name_col_1based: Optional[int] = None,
    strand_col_1based: Optional[int] = None,
) -> Dict[str, str]:
    """
    Load BED mapping name->strand. Supports:
      - chrom start end name score strand            (name=4, strand=6)
      - chrom start end strand name                  (name=5, strand=4)
    Or override with --bed-name-col / --bed-strand-col (1-based).
    """
    name_to_strand: Dict[str, str] = {}
    conflicts: List[Tuple[str, str, str]] = []

    forced = (name_col_1based is not None and strand_col_1based is not None)
    n_idx_forced = (name_col_1based - 1) if name_col_1based else None
    s_idx_forced = (strand_col_1based - 1) if strand_col_1based else None

    with smart_open(bed_path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            f = s.split()
            if len(f) < 5:
                raise ValueError(f"BED parse error on line {ln}: expected >=5 fields, got {len(f)}")

            # Determine columns
            if forced:
                try:
                    name = f[n_idx_forced]
                    strand = f[s_idx_forced]
                except IndexError:
                    raise ValueError(f"BED column index out of range on line {ln}")
            else:
                # Heuristics: accept two common layouts
                name = None
                strand = None
                if len(f) >= 6 and _is_strand_token(f[5]):
                    # chrom start end name score strand
                    name = f[3]
                    strand = f[5]
                elif len(f) >= 5 and _is_strand_token(f[3]):
                    # chrom start end strand name
                    name = f[4]
                    strand = f[3]
                else:
                    # Try another reasonable guess: if f[4] is '+'/'-': chrom start end ? strand ?
                    if len(f) >= 6 and _is_strand_token(f[4]):
                        # then assume name in f[5]
                        name = f[5]
                        strand = f[4]
                    else:
                        raise ValueError(
                            f"BED auto-detect failed on line {ln}. "
                            f"Provide --bed-name-col and --bed-strand-col."
                        )

            if not _is_strand_token(strand):
                raise ValueError(f"BED invalid strand '{strand}' on line {ln} (name={name})")

            prev = name_to_strand.get(name)
            if prev is None:
                name_to_strand[name] = strand
            elif prev != strand:
                conflicts.append((name, prev, strand))

    if conflicts:
        msg = "\n  ".join(f"name '{n}' has conflicting strands: '{a}' vs '{b}'" for n, a, b in conflicts)
        raise ValueError("BED contains conflicting strand assignments:\n  " + msg)

    if not name_to_strand:
        raise ValueError("BED contained no usable rows.")

    return name_to_strand

def parse_float(x: str, ctx: str) -> float:
    try:
        return float(x)
    except Exception as e:
        raise ValueError(f"Failed to parse float for {ctx}: {x!r}") from e

def stream_modkit_rows(path: str):
    # Modkit: col4=name, col5=mean entropy, col6=strand
    with smart_open(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            raw = line.rstrip("\n")
            if not raw or raw.startswith("#"):
                continue
            f = raw.split()
            if len(f) < 6:
                raise ValueError(f"Modkit parse error on line {ln} in {path}: expected >=6 fields, got {len(f)}")
            name = f[3]
            strand = f[5]
            if strand not in {"+", "-"}:
                raise ValueError(f"Modkit invalid strand '{strand}' on line {ln} in {path} (name={name})")
            mean_entropy = parse_float(f[4], f"mean_entropy (line {ln} in {path})")
            yield name, strand, mean_entropy

def main():
    ap = argparse.ArgumentParser(description="Summarize Modkit entropy across multiple inputs with strand enforcement via BED.")
    ap.add_argument("-i", "--input", action="append", required=True,
                    help="Modkit entropy file (TSV or .gz). Repeat for multiple files.")
    ap.add_argument("--bed", required=True, help="BED file (TSV or .gz).")
    ap.add_argument("--bed-name-col", type=int, help="1-based column index for NAME in BED (overrides auto-detect).")
    ap.add_argument("--bed-strand-col", type=int, help="1-based column index for STRAND in BED (overrides auto-detect).")
    ap.add_argument("--name", required=True, help="Label printed above the stats block.")
    ap.add_argument("--quiet", action="store_true", help="Suppress stderr summary.")
    args = ap.parse_args()

    try:
        name_to_strand = load_bed_name_to_strand(
            args.bed,
            name_col_1based=args.bed_name_col,
            strand_col_1based=args.bed_strand_col,
        )
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    included = 0
    excl_wrong_strand = 0
    excl_missing_name = 0
    entropies: List[float] = []
    files_processed = 0

    try:
        for ipath in args.input:
            files_processed += 1
            for name, row_strand, mean_entropy in stream_modkit_rows(ipath):
                bed_strand = name_to_strand.get(name)
                if bed_strand is None:
                    excl_missing_name += 1
                    continue
                if row_strand != bed_strand:
                    excl_wrong_strand += 1
                    continue
                included += 1
                entropies.append(mean_entropy)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    if included == 0:
        if not args.quiet:
            print("[ERROR] No rows matched the correct strand across provided inputs.", file=sys.stderr)
            print("        Verify Modkit col4=name, col6=strand and that names exist in BED.", file=sys.stderr)
        sys.exit(2)

    ent_min = min(entropies)
    ent_max = max(entropies)
    ent_mean = sum(entropies) / len(entropies)

    if not args.quiet:
        print("# Combined summary (after enforcing BED strand):", file=sys.stderr)
        print(f"  Files processed:           {files_processed}", file=sys.stderr)
        print(f"  Included rows:             {included}", file=sys.stderr)
        print(f"  Excluded (wrong strand):   {excl_wrong_strand}", file=sys.stderr)
        print(f"  Excluded (name not in BED):{excl_missing_name}", file=sys.stderr)
        print(f"  Entropy range (min..max):  {ent_min:.9f} .. {ent_max:.9f}", file=sys.stderr)
        print(f"  Mean of means (entropy):   {ent_mean:.9f}", file=sys.stderr)

    # Append-friendly 3-line block
    print(args.name)
    print("\t".join([
        "included_rows",
        "excluded_wrong_strand",
        "excluded_missing_name",
        "files_processed",
        "entropy_min",
        "entropy_max",
        "entropy_mean",
    ]))
    print("\t".join([
        str(included),
        str(excl_wrong_strand),
        str(excl_missing_name),
        str(files_processed),
        f"{ent_min:.9f}",
        f"{ent_max:.9f}",
        f"{ent_mean:.9f}",
    ]))

if __name__ == "__main__":
    main()
