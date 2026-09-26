#!/usr/bin/env python3
import argparse, sys, gzip, re
from pathlib import Path
from collections import defaultdict

def open_text(p: Path):
    return gzip.open(p, "rt", errors="ignore") if str(p).endswith(".gz") else open(p, "rt", errors="ignore")

def load_interest(path: Path):
    keep = set()
    with open(path, "rt", errors="ignore") as fh:
        for i, ln in enumerate(fh):
            if not ln.strip(): 
                continue
            tok = ln.strip().split()
            if not tok: 
                continue
            # skip header
            if i == 0 and (tok[0].lower() in ("locus", "gene", "locus_tag")):
                continue
            keep.add(tok[0])
    return keep

def load_bed(path: Path):
    bed = defaultdict(list)  # contig -> list of (start,end,locus)
    with open(path, "rt", errors="ignore") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"): 
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            locus = parts[3] if len(parts) >= 4 and parts[3] else f"{chrom}:{start}-{end}"
            bed[chrom].append((start, end, locus))
    for chrom in bed:
        bed[chrom].sort()
    return bed

def variant_type(ref: str, alt_field: str):
    alts = alt_field.split(",")
    if len(ref) == 1 and all(len(a) == 1 for a in alts):
        return "snp"
    return "indel"

def alt_max_len(ref: str, alt_field: str):
    alts = alt_field.split(",")
    return max(abs(len(a) - len(ref)) for a in alts)

def parse_vcf_records(vcf_path: Path, vtype: str, min_qual: float):
    with open_text(vcf_path) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, pos, _id, ref, alt, qual = parts[:6]
            if qual != ".":
                try:
                    if float(qual) < min_qual:
                        continue
                except Exception:
                    pass
            vclass = variant_type(ref, alt)
            if vtype == "snp" and vclass != "snp":
                continue
            if vtype == "indel" and vclass != "indel":
                continue
            if vtype == "long_indel":
                if vclass != "indel":
                    continue
                if alt_max_len(ref, alt) <= 50:
                    continue
            yield chrom, int(pos), ref, alt

def count_overlaps_for_sample(vcf_path: Path, bed, vtype: str, min_qual: float):
    counts = defaultdict(int)
    for chrom, pos, ref, alt in parse_vcf_records(vcf_path, vtype=vtype, min_qual=min_qual):
        if chrom not in bed:
            continue
        for start, end, locus in bed[chrom]:
            if pos < start:
                break
            if start <= pos <= end:
                counts[locus] += 1
    return counts

def infer_group(sample: str):
    m = re.match(r"^([A-Za-z]+)(T|C)\d+$", sample)
    if m:
        letters, tc = m.groups()
        return letters.upper(), ("T" if tc.upper() == "T" else "C")
    return "ALL", "T"

def load_group_map(path: Path):
    mp = {}
    with open(path, "rt", errors="ignore") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"): 
                continue
            toks = ln.rstrip("\n").split("\t")
            if len(toks) < 2:
                continue
            sample = toks[0]
            group = toks[1]
            status = toks[2].upper() if len(toks) >= 3 else "T"
            mp[sample] = (group, status)
    return mp

def discover_vcfs(vcf_dir: Path, vcf_name: str):
    mapping = {}
    for sample_dir in sorted(p for p in vcf_dir.iterdir() if p.is_dir()):
        sample = sample_dir.name
        vcf = sample_dir / vcf_name
        if vcf.exists():
            mapping[sample] = vcf
    return mapping

def main():
    ap = argparse.ArgumentParser(description="Locus x treatment matrix of variant counts from Clair3 outputs.")
    ap.add_argument("--bed", required=True, help="Gene BED (chrom start end locusTag ...)")
    ap.add_argument("--in", dest="interest", required=True, help="File with locus tags of interest (first column).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--vcf-dir", help="Root folder with per-sample Clair3 outputs (subfolders = samples).")
    g.add_argument("--vcf", nargs="+", help="One or more VCF paths (per-sample).")
    ap.add_argument("--vcf-name", default="merge_output.vcf.gz", help="VCF name under each sample dir (default: merge_output.vcf.gz)")
    ap.add_argument("--group-map", help="Optional TSV: sample<TAB>group<TAB>status(T/C). If omitted, infer from sample name.")
    ap.add_argument("--type", choices=["all","snp","indel","long_indel"], default="all", help="Variant type to count (default: all)")
    ap.add_argument("--min-qual", type=float, default=0.0, help="Minimum site QUAL to include (default: 0)")
    ap.add_argument("--treat-only", action="store_true", help="Only aggregate across treatment samples (status T)")
    args = ap.parse_args()

    bed = load_bed(Path(args.bed))
    interest = load_interest(Path(args.interest))

    # Gather VCFs
    sample_to_vcf = {}
    if args.vcf_dir:
        sample_to_vcf = discover_vcfs(Path(args.vcf_dir), args.vcf_name)
    else:
        for pstr in args.vcf:
            p = Path(pstr)
            if p.is_dir():
                v = p / args.vcf_name
                if v.exists():
                    sample = p.name
                    sample_to_vcf[sample] = v
            else:
                sample = p.parent.name
                sample_to_vcf[sample] = p

    if not sample_to_vcf:
        print("No VCFs found; check --vcf-dir or --vcf inputs.", file=sys.stderr)
        sys.exit(1)

    # Group map and columns
    group_map = {}
    if args.group_map:
        group_map = load_group_map(Path(args.group_map))

    groups = []
    sample_meta = {}
    for sample in sample_to_vcf:
        if sample in group_map:
            grp, status = group_map[sample]
        else:
            grp, status = infer_group(sample)
        sample_meta[sample] = (grp, status)
        if grp not in groups:
            groups.append(grp)

    # Init matrix
    matrix = {locus: {grp: 0 for grp in groups} for locus in interest}

    # Accumulate counts
    for sample, vcf in sample_to_vcf.items():
        grp, status = sample_meta[sample]
        if args.treat_only and status != "T":
            continue
        per_locus = count_overlaps_for_sample(Path(vcf), bed, vtype=args.type, min_qual=args.min_qual)
        for locus, n in per_locus.items():
            if locus in matrix:
                matrix[locus][grp] += n

    # Emit
    header = ["locus"] + groups
    print("\t".join(header))
    for locus in sorted(matrix.keys()):
        row = [locus] + [str(matrix[locus].get(grp, 0)) for grp in groups]
        print("\t".join(row))

if __name__ == "__main__":
    main()
