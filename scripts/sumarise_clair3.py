#!/usr/bin/env python3
"""
Summarise Clair3 outputs per sample *per contig*.

- Input root (--in): directory containing per-sample Clair3 output folders.
  Each sample folder should contain a VCF (default: merge_output.vcf.gz).
- FAI (-fai): reference .fai file (contig lengths) to compute variant densities.
- Output: TSV to stdout with columns:
    sample  contig  contig_len  n_variants  n_snps  n_indels  mean_AF  median_QUAL  var_per_kb  var_per_mb

Notes:
- AF is taken from INFO/AF if present; otherwise left blank.
- Variant type: SNP if len(REF)==1 and all ALT alleles have length 1; otherwise INDEL.
- Works with .vcf or .vcf.gz.
"""

import sys, re, gzip
from pathlib import Path
import argparse
import statistics as stats

def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")

def parse_fai(fai_path: Path):
    lengths = {}
    if fai_path and fai_path.exists():
        with open(fai_path, "rt", errors="ignore") as fh:
            for ln in fh:
                if not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
    return lengths

def is_snp(ref: str, alt_field: str) -> bool:
    # alt_field may contain comma-separated ALTs
    alts = alt_field.split(",")
    if len(ref) != 1:
        return False
    for a in alts:
        if len(a) != 1:
            return False
    return True

def parse_info_af(info: str):
    # Return list of AF floats if present, else empty list
    # Look for AF=0.12 or AF=0.1,0.2
    for field in info.split(";"):
        if field.startswith("AF="):
            try:
                vals = [float(x) for x in field[3:].split(",") if x not in (".","")]
                return vals
            except Exception:
                return []
    return []

def summarise_vcf(vcf_path: Path, contig_lengths: dict):
    # returns dict: contig -> stats
    # stats: {"n_vars":int,"n_snps":int,"n_indels":int,"AFs":[float...],"QUALs":[float...]}
    out = {}
    if not vcf_path.exists():
        return out
    with open_text(vcf_path) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, pos, _id, ref, alt, qual = parts[:6]
            info = parts[7] if len(parts) > 7 else ""
            # Build entry
            st = out.setdefault(chrom, {"n_vars":0,"n_snps":0,"n_indels":0,"AFs":[], "QUALs":[]})
            st["n_vars"] += 1
            if is_snp(ref, alt):
                st["n_snps"] += 1
            else:
                st["n_indels"] += 1
            # QUAL may be '.'
            try:
                if qual != ".":
                    st["QUALs"].append(float(qual))
            except Exception:
                pass
            AFs = parse_info_af(info)
            st["AFs"].extend(AFs)
    return out

def main():
    ap = argparse.ArgumentParser(description="Per-contig Clair3 summary to TSV (stdout).")
    ap.add_argument("--in", dest="root", required=True, help="Root folder containing per-sample Clair3 outputs (subfolders = samples).")
    ap.add_argument("--fai", dest="fai", required=True, help="Reference .fai file for contig lengths.")
    ap.add_argument("--vcf-name", default="merge_output.vcf.gz", help="VCF file name (default: merge_output.vcf.gz).")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    fai = Path(args.fai).resolve()
    contig_lengths = parse_fai(fai)

    # Header
    cols = ["sample","contig","contig_len","n_variants","n_snps","n_indels","mean_AF","median_QUAL","var_per_kb","var_per_mb"]
    print("\t".join(cols))

    # Iterate subdirs
    for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        sample = sample_dir.name
        vcf = sample_dir / args.vcf_name
        stats_by_contig = summarise_vcf(vcf, contig_lengths)

        # Emit rows per contig
        for ctg, st in sorted(stats_by_contig.items()):
            clen = contig_lengths.get(ctg, 0)
            n = st["n_vars"]
            per_kb = (n / clen * 1e3) if clen > 0 else ""
            per_mb = (n / clen * 1e6) if clen > 0 else ""
            mean_af = stats.mean(st["AFs"]) if st["AFs"] else ""
            median_qual = stats.median(st["QUALs"]) if st["QUALs"] else ""
            row = [
                sample,
                ctg,
                str(clen if clen else ""),
                str(n),
                str(st["n_snps"]),
                str(st["n_indels"]),
                f"{mean_af:.4f}" if mean_af != "" else "",
                f"{median_qual:.2f}" if median_qual != "" else "",
                f"{per_kb:.6f}" if per_kb != "" else "",
                f"{per_mb:.6f}" if per_mb != "" else "",
            ]
            print("\t".join(row))

        # Add a TOTAL row per sample
        total_n = sum(st["n_vars"] for st in stats_by_contig.values())
        total_snps = sum(st["n_snps"] for st in stats_by_contig.values())
        total_indels = sum(st["n_indels"] for st in stats_by_contig.values())
        # weighted mean AF is ambiguous; report simple mean over all AFs
        all_afs = [af for st in stats_by_contig.values() for af in st["AFs"]]
        all_quals = [q for st in stats_by_contig.values() for q in st["QUALs"]]
        mean_af_total = stats.mean(all_afs) if all_afs else ""
        median_qual_total = stats.median(all_quals) if all_quals else ""
        row = [
            sample,
            "__TOTAL__",
            "",
            str(total_n),
            str(total_snps),
            str(total_indels),
            f"{mean_af_total:.4f}" if mean_af_total != "" else "",
            f"{median_qual_total:.2f}" if median_qual_total != "" else "",
            "",
            "",
        ]
        print("\t".join(row))

if __name__ == "__main__":
    main()
