#!/usr/bin/env python3
"""
entropy_stats_calc.py

Terminal-only summary statistics for methylation entropy changes.
Preserves the full input shell and gene/operon logic of entropy_express_calc.py,
but replaces all plotting and regression with concise reporting.

Reports, per treatment and per strand context:
- mean |Δ-entropy|
- standard deviation and variance
- number and fraction of ~zero-change sites
- minimum, maximum, and largest |Δ-entropy|
"""

import argparse
import sys
import numpy as np
import pandas as pd
from collections import Counter

ZERO_EPS = 0.005    # threshold for "approximately zero" entropy change

# ─────────── I/O helpers ───────────
def read_bed(path):
    df = pd.read_csv(
        path, sep="\t", header=None, usecols=[3, 4, 5],
        names=["locus_tag", "entropy", "strand"]
    )
    df["locus_tag"] = df["locus_tag"].astype(str)
    df["strand"] = df["strand"].astype(str)
    df["entropy"] = pd.to_numeric(df["entropy"], errors="coerce")
    df = df[df["strand"].isin(["+", "-"])].dropna(subset=["entropy"])
    return df.reset_index(drop=True)

def read_de(path):
    df = pd.read_csv(path, sep=None, engine="python")
    df.rename(columns={df.columns[0]: "locus_tag"}, inplace=True)
    lfc_col = None
    for c in df.columns:
        if c.lower() in ("log2foldchange", "log2fc", "log2_fc", "log2_fold_change"):
            lfc_col = c
            break
    if lfc_col is None:
        raise SystemExit(f"{path}: no log2FC column found.")
    return df.set_index("locus_tag")[lfc_col].astype(float).to_dict()

def read_gff(path):
    strands = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            strand = parts[6]
            if strand not in {"+", "-"}:
                continue
            for attr in parts[8].split(";"):
                if attr.startswith("locus_tag="):
                    strands[attr.split("=", 1)[1]] = strand
                    break
    if not strands:
        raise SystemExit("No locus_tag entries parsed from GFF.")
    return strands

def read_operon_map(path):
    m = pd.read_csv(path, sep="\t", header=None, names=["operon_id", "locus_tag"])
    m["operon_id"] = m["operon_id"].astype(str)
    m["locus_tag"] = m["locus_tag"].astype(str)
    return m.groupby("operon_id")["locus_tag"].apply(list).to_dict()

# ─────────── mode detection ───────────
def detect_mode(T_bed, C_bed, gff_strand):
    gene_tags = set(gff_strand.keys())
    all_tags = set()
    for d in (T_bed, C_bed):
        for df in d.values():
            all_tags.update(df["locus_tag"].unique())
    frac = sum(t in gene_tags for t in all_tags) / max(1, len(all_tags))
    return "gene" if frac >= 0.5 else "operon"

def operon_strand_map(op2genes, gff_strand):
    op2strand = {}
    for op, genes in op2genes.items():
        strands = [gff_strand[g] for g in genes if g in gff_strand]
        if strands:
            op2strand[op] = Counter(strands).most_common(1)[0][0]
    return op2strand

# ─────────── table builders ───────────
def build_gene_tables(gff_strand, T_bed, C_bed, trts):
    rows_m, rows_o = [], []
    for tag, strand in gff_strand.items():
        opp = "-" if strand == "+" else "+"
        rm, ro = {"locus_tag": tag}, {"locus_tag": tag}
        for t in trts:
            tm = T_bed[t].query(
                "locus_tag == @tag and strand == @strand"
            )["entropy"].mean()
            cm = C_bed[t].query(
                "locus_tag == @tag and strand == @strand"
            )["entropy"].mean()
            to = T_bed[t].query(
                "locus_tag == @tag and strand == @opp"
            )["entropy"].mean()
            co = C_bed[t].query(
                "locus_tag == @tag and strand == @opp"
            )["entropy"].mean()

            rm[f"Delta_{t}"] = abs(tm - cm) if pd.notna(tm) and pd.notna(cm) else np.nan
            ro[f"Delta_{t}"] = abs(to - co) if pd.notna(to) and pd.notna(co) else np.nan

        rows_m.append(rm)
        rows_o.append(ro)

    return pd.DataFrame(rows_m), pd.DataFrame(rows_o)

def build_operon_tables(T_bed, C_bed, trts, op2genes, gff_strand):
    """
    Operon entropy is read DIRECTLY from BED files (which use operon IDs).
    The operon→gene map is used ONLY to infer coding strand via GFF.
    """
    # infer coding strand per operon from member genes
    op2strand = {}
    for op, genes in op2genes.items():
        strands = [gff_strand[g] for g in genes if g in gff_strand]
        if strands:
            op2strand[op] = Counter(strands).most_common(1)[0][0]

    rows_m, rows_o = [], []

    for op, strand in op2strand.items():
        opp = "-" if strand == "+" else "+"
        rm = {"locus_tag": f"operon_{op}"}
        ro = {"locus_tag": f"operon_{op}"}

        for t in trts:
            tm = T_bed[t].query(
                "locus_tag == @op and strand == @strand"
            )["entropy"].mean()
            cm = C_bed[t].query(
                "locus_tag == @op and strand == @strand"
            )["entropy"].mean()

            to = T_bed[t].query(
                "locus_tag == @op and strand == @opp"
            )["entropy"].mean()
            co = C_bed[t].query(
                "locus_tag == @op and strand == @opp"
            )["entropy"].mean()

            rm[f"Delta_{t}"] = abs(tm - cm) if pd.notna(tm) and pd.notna(cm) else np.nan
            ro[f"Delta_{t}"] = abs(to - co) if pd.notna(to) and pd.notna(co) else np.nan

        rows_m.append(rm)
        rows_o.append(ro)

    return pd.DataFrame(rows_m), pd.DataFrame(rows_o)


# ─────────── reporting ───────────
def report(df, label):
    print(f"\n=== {label} strand ===")

    delta_cols = sorted(c for c in df.columns if c.startswith("Delta_"))
    if not delta_cols:
        print("No Δ-entropy columns present.")
        return

    for col in delta_cols:
        t = col.replace("Delta_", "")
        x = df[col].dropna().values
        if len(x) == 0:
            print(f"\nTreatment: {t}")
            print("  no valid entropy deltas")
            continue

        nz = np.sum(x < ZERO_EPS)

        print(f"\nTreatment: {t}")
        print(f"  n sites              : {len(x)}")
        print(f"  mean |Δ-entropy|     : {np.mean(x):.6f}")
        print(f"  SD                   : {np.std(x, ddof=1):.6f}")
        print(f"  variance             : {np.var(x, ddof=1):.6f}")
        print(f"  min / max            : {np.min(x):.6f} / {np.max(x):.6f}")
        print(f"  max |Δ|              : {np.max(x):.6f}")
        print(f"  ~zero-change sites   : {nz} ({100*nz/len(x):.2f}%)")

# ─────────── main ───────────
def main():
    pa = argparse.ArgumentParser()
    pa.add_argument("-i", action="append", nargs=2, required=True)
    pa.add_argument("-c", action="append", nargs=2, required=True)
    pa.add_argument("-d", action="append", nargs=2, required=True)
    pa.add_argument("--operon", required=False)
    pa.add_argument("--gff", required=True)
    pa.add_argument("-o", "--outdir", required=True)
    pa.add_argument("-p", "--prefix", required=True)
    args = pa.parse_args()

    gff_strand = read_gff(args.gff)
    T_bed = {n: read_bed(p) for p, n in args.i}
    C_bed = {n: read_bed(p) for p, n in args.c}
    trts = sorted(T_bed.keys())

    mode = detect_mode(T_bed, C_bed, gff_strand)
    print(f"[INFO] Analysis mode: {mode}")

    if mode == "gene":
        df_m, df_o = build_gene_tables(gff_strand, T_bed, C_bed, trts)
    else:
        if not args.operon:
            raise SystemExit("Operon mode detected but --operon not provided.")
        op2genes = read_operon_map(args.operon)
        df_m, df_o = build_operon_tables(T_bed, C_bed, trts, op2genes, gff_strand)

    print("\n==============================")
    print("Entropy change summary report")
    print("==============================")
    report(df_m, "Matching")
    report(df_o, "Opposing")

if __name__ == "__main__":
    main()
