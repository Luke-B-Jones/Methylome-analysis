#!/usr/bin/env python3
"""
entropy_express_calc.py — |Δ-entropy| vs |log2FC| (matching/opposing strands), with ordered percentile bins,
bootstrap CIs for OLS slopes, and permutation-null histograms.

Key behaviors:
- Bins are defined by PERCENTILES and labeled by order (Bin 1, Bin 2, …), not hard-coded percentiles.
- Δ = treatment - control (for each treatment name).
- Operon mapping (2-col TSV: operon_id<TAB>locus_tag) is used ONLY if BED locus_tags
  are not recognized gene tags from the GFF (e.g., numeric operon IDs).
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from collections import Counter
from matplotlib.ticker import MaxNLocator

# ─────────── User-configurable params ───────────
PERCENTILES      = [50, 75, 95]   # cutoff percentiles (define len+1 bins)
BOOTSTRAP_B      = 500            # bootstrap replicates
BOOTSTRAP_ALPHA  = 0.05           # 95% CI

# A simple palette; first bin is grey/black, the rest cycle through this palette:
PALETTE          = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00", "#56B4E9"]
FIRST_BIN_COLOR  = "black"
FIRST_BIN_ALPHA  = 0.10
BIN_ALPHA_HIGH   = 0.45           # alpha for higher bins
FIRST_BIN_SIZE   = 4
OTHER_BIN_SIZE   = 7

# fonts (global)
plt.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 24,
    "axes.labelsize": 24,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})

# ─────────── Helpers: I/O ───────────
def read_bed(path):
    """
    Expect a BED-like file with at least 6 columns; use columns 4,5,6 (0-based: 3,4,5):
      col4 = locus_tag (string: gene tag or operon ID),
      col5 = entropy (numeric),
      col6 = strand ('+'/'-').
    """
    try:
        df = pd.read_csv(path, sep="\t", header=None, usecols=[3,4,5],
                         names=["locus_tag","entropy","strand"])
    except Exception as e:
        raise SystemExit(f"BED parse error for {path}: {e}")
    df["locus_tag"] = df["locus_tag"].astype(str)
    df["strand"]    = df["strand"].astype(str)
    df["entropy"]   = pd.to_numeric(df["entropy"], errors="coerce")
    df = df[df["strand"].isin(["+","-"])].dropna(subset=["entropy"]).reset_index(drop=True)
    return df

def read_de(path):
    """
    CSV/TSV with first column locus_tag and a log2FC column.
    Recognized (case-insensitive): 'log2FoldChange','log2fc','log2_fc','log2_fold_change'.
    Returns dict: locus_tag -> float(log2FC).
    """
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as e:
        raise SystemExit(f"DE file parse error for {path}: {e}")
    if df.shape[1] < 2:
        raise SystemExit(f"DE CSV {path} must have at least two columns (locus_tag, log2FC).")
    df.rename(columns={df.columns[0]: "locus_tag"}, inplace=True)
    lfc_col = None
    for c in df.columns:
        cl = c.lower()
        if cl in ("log2foldchange","log2fc","log2_fc","log2_fold_change"):
            lfc_col = c
            break
    if lfc_col is None:
        raise SystemExit(
            f"DE file {path} missing a log2FC column (accepted: log2FoldChange/log2FC/log2_fc/log2_fold_change)."
        )
    return df.set_index("locus_tag")[lfc_col].astype(float).to_dict()

def read_gff(path):
    """Map locus_tag -> strand from GFF. Requires 'locus_tag=' in attributes."""
    strands = {}
    with open(path) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            strand = parts[6]
            if strand not in {"+","-"}:
                continue
            for attr in parts[8].split(";"):
                if attr.startswith("locus_tag="):
                    tag = attr.split("=",1)[1]
                    strands[tag] = strand
                    break
    if not strands:
        raise SystemExit(f"No locus_tag entries parsed from GFF: {path}")
    return strands

def read_operon_map(path):
    """
    2-column TSV: operon_id<TAB>locus_tag (one row per gene).
    Returns:
      op2genes: dict operon_id -> list of locus_tags
    """
    try:
        m = pd.read_csv(path, sep="\t", header=None, names=["operon_id","locus_tag"], usecols=[0,1])
    except Exception as e:
        raise SystemExit(f"Operon map parse error for {path}: {e}")
    m["operon_id"]  = m["operon_id"].astype(str)
    m["locus_tag"]  = m["locus_tag"].astype(str)
    op2genes = m.groupby("operon_id")["locus_tag"].apply(list).to_dict()
    if not op2genes:
        raise SystemExit(f"Operon map is empty: {path}")
    return op2genes

# ─────────── Helpers: logic ───────────
def gen_bin_styles(n_bins):
    """Return colors, alphas, sizes arrays of length n_bins."""
    colors = [FIRST_BIN_COLOR] + [PALETTE[i % len(PALETTE)] for i in range(n_bins - 1)]
    alphas = [FIRST_BIN_ALPHA] + [BIN_ALPHA_HIGH] * (n_bins - 1)
    sizes  = [FIRST_BIN_SIZE]  + [OTHER_BIN_SIZE] * (n_bins - 1)
    return colors, alphas, sizes

def has_variation(x):
    return (x.size >= 3) and (np.nanvar(x) > 0)

def bootstrap_ci(x, y, B=BOOTSTRAP_B, alpha=BOOTSTRAP_ALPHA):
    n = len(x)
    betas, r2s = [], []
    for _ in range(B):
        idx = np.random.randint(0, n, n)
        xb, yb = x[idx], y[idx]
        if not has_variation(xb):
            continue
        m = sm.OLS(yb, sm.add_constant(xb)).fit()
        betas.append(float(m.params[1]))
        r2s.append(float(m.rsquared))
    if len(betas) == 0:
        return (np.array([np.nan, np.nan]), np.array([np.nan, np.nan]))
    lo, hi = 100*alpha/2, 100*(1 - alpha/2)
    return np.percentile(betas, [lo, hi]), np.percentile(r2s, [lo, hi])

def bin_edges_from_percentiles(values, percentiles):
    thr = np.percentile(values, sorted(percentiles))
    return [0.0] + thr.tolist() + [np.inf]

def select_bin_mask(x, edges, i):
    # Bin i in 0..n-1: i=0 => [lo,hi], i>0 => (lo,hi]
    lo, hi = edges[i], edges[i+1]
    return ((x >= lo) & (x <= hi)) if i == 0 else ((x > lo) & (x <= hi))

def detect_mode(t_beds, c_beds, gff_strand):
    """
    Decide 'gene' vs 'operon' mode by checking if BED locus_tags match GFF locus_tags.
    If <50% of unique tags across all BEDs are in GFF, assume operon mode.
    """
    gene_tags = set(gff_strand.keys())
    all_tags = set()
    for d in (t_beds, c_beds):
        for df in d.values():
            all_tags.update(df["locus_tag"].unique().tolist())
    if not all_tags:
        raise SystemExit("No locus_tag entries found in BED inputs.")
    in_gff = sum(1 for t in all_tags if t in gene_tags)
    frac = in_gff / max(1, len(all_tags))
    return "gene" if frac >= 0.5 else "operon"

def operon_strand_map(op2genes, gff_strand):
    """
    For each operon_id, infer strand from member genes.
    If mixed, use majority and print a warning.
    """
    op2strand = {}
    for op, genes in op2genes.items():
        strands = [gff_strand.get(g) for g in genes if g in gff_strand]
        strands = [s for s in strands if s in {"+","-"}]
        if not strands:
            continue
        cnt = Counter(strands)
        if len(cnt) > 1:
            # choose majority strand
            sys.stderr.write(f"[WARN] Operon {op} has mixed strands {dict(cnt)}; using majority.\n")
        op2strand[op] = cnt.most_common(1)[0][0]
    return op2strand

# ─────────── analysis & plotting ───────────
def scatter_percentiles_with_ci(df, trts, stem, outdir):
    # Build long-form pooled across treatments
    frames = []
    for t in trts:
        a = df.get(f"Delta_{t}", pd.Series(dtype=float)).abs()
        b = df.get(f"Log2FC_{t}", pd.Series(dtype=float)).abs()
        tmp = pd.DataFrame({"abs_delta": a, "abs_lfc": b}).dropna()
        if not tmp.empty:
            frames.append(tmp)
    if not frames:
        print(f"[WARN] No data for {stem}; skipping.")
        return
    dat = pd.concat(frames, ignore_index=True)

    # Output folder
    figdir = os.path.join(outdir, "percentile_analysis")
    os.makedirs(figdir, exist_ok=True)

    # Edges and dynamic styles
    edges = bin_edges_from_percentiles(dat.abs_delta.values, PERCENTILES)
    n_bins = len(edges) - 1
    COLORS, ALPHAS, SIZES = gen_bin_styles(n_bins)

    # ── Scatter by bins ──
    plt.figure(figsize=(8, 6))
    xvals = dat.abs_delta.values
    for i in range(n_bins):
        sub = dat[select_bin_mask(xvals, edges, i)]
        if sub.empty:
            continue
        plt.scatter(sub.abs_delta, sub.abs_lfc,
                    s=SIZES[i], c=COLORS[i], alpha=ALPHAS[i], linewidths=0)

    plt.xlabel(r"$|z(\Delta\mathrm{ME})|$")
    plt.ylabel(r"$|\mathrm{Log}_2\mathrm{FC}|$")

    # Axis scaling
    x_max = float(dat.abs_delta.max())
    x_tick_max = np.floor(x_max * 100) / 100.0
    x_plot_max = max(0.05, x_tick_max + 0.01)
    xticks = np.round(np.linspace(0, x_tick_max, 4), 3)

    y_max_data = float(dat.abs_lfc.max())
    y_step = 3.0
    y_max = y_step * np.ceil(max(6.0, y_max_data) / y_step)
    yticks = np.arange(0.0, y_max + 1e-9, y_step)

    plt.xlim(0, x_plot_max)
    plt.ylim(0, y_max)
    plt.xticks(xticks)
    plt.yticks(yticks)

    ax = plt.gca()
    ax.set_xticklabels(['0'] + [f"{x:.3f}" for x in xticks[1:]])
    ax.set_yticklabels(['0'] + [str(int(y)) for y in yticks[1:]])

    # ── Fits & bootstrap CIs for higher bins (bins 2..n_bins) ──
    xfit = np.linspace(0, x_plot_max, 300)
    for j in range(1, n_bins):  # skip bin 0
        sub = dat[select_bin_mask(xvals, edges, j)]
        if sub.shape[0] < 3 or not has_variation(sub.abs_delta.values):
            continue
        x = sub.abs_delta.values
        y = sub.abs_lfc.values
        model = sm.OLS(y, sm.add_constant(x)).fit()
        beta, alpha0 = float(model.params[1]), float(model.params[0])
        pval, R2 = float(model.pvalues[1]), float(model.rsquared)
        (beta_lo, beta_hi), (r2_lo, r2_hi) = bootstrap_ci(x, y)

        # CI band & line
        y_lo = alpha0 + beta_lo * xfit if np.isfinite(beta_lo) else np.full_like(xfit, np.nan)
        y_hi = alpha0 + beta_hi * xfit if np.isfinite(beta_hi) else np.full_like(xfit, np.nan)
        if np.all(np.isfinite(y_lo)) and np.all(np.isfinite(y_hi)):
            plt.fill_between(xfit, y_lo, y_hi, color=COLORS[j], alpha=0.2)
        plt.plot(xfit, alpha0 + beta * xfit, color=COLORS[j], lw=2)

        # Annotation (by bin order)
        n = len(x)
        pct = PERCENTILES[j-1]
        txt = (
            f"{pct}th (n={n})  p={pval:.1e}\n"
            f"β={beta:.2f} [{np.nan_to_num(beta_lo):.2f},{np.nan_to_num(beta_hi):.2f}]  "
            f"R²={R2:.2f} [{np.nan_to_num(r2_lo):.2f},{np.nan_to_num(r2_hi):.2f}]"
        )

        top_y = 0.98
        line_spacing = 0.15
        plt.text(0.97, top_y - (j-1)*line_spacing, txt,
                 ha="right", va="top", transform=ax.transAxes, color=COLORS[j])

    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"{stem}_percentile_scatter.png"), dpi=300)
    plt.close()

    # ────────── Permutation-null histograms for higher bins ──────────
    rng = np.random.default_rng(1)
    plt.figure(figsize=(6, 4))
    for j in range(1, n_bins):
        sub = dat[select_bin_mask(xvals, edges, j)]
        if sub.shape[0] < 3 or not has_variation(sub.abs_delta.values):
            continue
        obs_beta = sm.OLS(sub.abs_lfc, sm.add_constant(sub.abs_delta)).fit().params.iloc[1]
        N_PERM = 1000
        null_b = np.empty(N_PERM)
        for k in range(N_PERM):
            x_shuff = rng.permutation(sub.abs_delta.values)
            null_b[k] = sm.OLS(sub.abs_lfc, sm.add_constant(x_shuff)).fit().params.iloc[1]
        plt.hist(null_b, bins=40, density=True, color=COLORS[j], alpha=0.30,
                 histtype='stepfilled', edgecolor=None)
        pct = PERCENTILES[j-1]
        plt.axvline(float(obs_beta), color=COLORS[j], lw=2, label=f"{pct}th")

    plt.xlabel("β (slope)")
    plt.ylabel("Density")
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
    handles, labels = ax.get_legend_handles_labels()
    legend_png = os.path.join(figdir, f"{stem}_slope_null_hist_legend.png")
    # Save legend separately
    fig_leg = plt.figure(figsize=(0.8 * max(1, len(labels)), 0.6))
    fig_leg.legend(handles, labels, loc="center", frameon=False, ncol=min(3, len(labels)),
                   fontsize=plt.rcParams["legend.fontsize"])
    fig_leg.savefig(legend_png, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig_leg)

    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"{stem}_slope_null_hist.png"), dpi=300)
    plt.close()

# ─────────── Main data assembly ───────────
def build_gene_tables(gff_strand, T_bed, C_bed, de_map, trts):
    rows_m, rows_o = [], []
    for tag, strand in gff_strand.items():
        opp = "-" if strand == "+" else "+"
        rm = {"locus_tag": tag}
        ro = {"locus_tag": tag}
        for n in trts:
            tdf, cdf = T_bed[n], C_bed[n]
            tm = tdf[(tdf.locus_tag == tag) & (tdf.strand == strand)].entropy.mean()
            cm = cdf[(cdf.locus_tag == tag) & (cdf.strand == strand)].entropy.mean()
            to = tdf[(tdf.locus_tag == tag) & (tdf.strand == opp)].entropy.mean()
            co = cdf[(cdf.locus_tag == tag) & (cdf.strand == opp)].entropy.mean()
            rm[f"Delta_{n}"]  = np.nan if np.isnan(tm) or np.isnan(cm) else float(tm - cm)   # Δ = T - C
            ro[f"Delta_{n}"]  = np.nan if np.isnan(to) or np.isnan(co) else float(to - co)   # Δ = T - C
            l2 = de_map[n].get(tag, np.nan)
            rm[f"Log2FC_{n}"] = float(l2) if pd.notna(l2) else np.nan
            ro[f"Log2FC_{n}"] = float(l2) if pd.notna(l2) else np.nan
        rows_m.append(rm); rows_o.append(ro)
    return (pd.DataFrame(rows_m).reset_index(drop=True),
            pd.DataFrame(rows_o).reset_index(drop=True))

def build_operon_tables(T_bed, C_bed, de_map, trts, op2genes, gff_strand):
    # infer per-operon strand (majority of member genes)
    op2strand = operon_strand_map(op2genes, gff_strand)
    # collect all operon IDs present in either T or C inputs
    all_ops = set()
    for d in (T_bed, C_bed):
        for df in d.values():
            all_ops.update(df["locus_tag"].unique().tolist())
    rows_m, rows_o = [], []
    for op in sorted(all_ops):
        op_strand = op2strand.get(op)
        if op_strand not in {"+","-"}:
            # cannot determine match/opposing strand for this operon (no gene strands); skip
            continue
        opp = "-" if op_strand == "+" else "+"
        rm = {"locus_tag": f"operon_{op}"}
        ro = {"locus_tag": f"operon_{op}"}
        for n in trts:
            tdf, cdf = T_bed[n], C_bed[n]
            tm = tdf[(tdf.locus_tag == op) & (tdf.strand == op_strand)].entropy.mean()
            cm = cdf[(cdf.locus_tag == op) & (cdf.strand == op_strand)].entropy.mean()
            to = tdf[(tdf.locus_tag == op) & (tdf.strand == opp)].entropy.mean()
            co = cdf[(cdf.locus_tag == op) & (cdf.strand == opp)].entropy.mean()
            rm[f"Delta_{n}"]  = np.nan if np.isnan(tm) or np.isnan(cm) else float(tm - cm)   # Δ = T - C
            ro[f"Delta_{n}"]  = np.nan if np.isnan(to) or np.isnan(co) else float(to - co)   # Δ = T - C
            # average log2FC across member genes
            genes = op2genes.get(op, [])
            vals = [de_map[n].get(g) for g in genes if g in de_map[n]]
            l2 = float(np.nanmean(vals)) if len(vals) else np.nan
            rm[f"Log2FC_{n}"] = l2
            ro[f"Log2FC_{n}"] = l2
        rows_m.append(rm); rows_o.append(ro)
    return (pd.DataFrame(rows_m).reset_index(drop=True),
            pd.DataFrame(rows_o).reset_index(drop=True))

def main():
    pa = argparse.ArgumentParser(
        description="Relate |Δ-entropy| (Δ = treatment - control) to |log2FC| across percentile bins (match/opposing strands), with CIs and permutation nulls."
    )
    pa.add_argument("-i", action="append", nargs=2, metavar=("BED","NAME"),
                    required=True, help="Treatment BED + name (Δ uses this as treatment)")
    pa.add_argument("-c", action="append", nargs=2, metavar=("BED","NAME"),
                    required=True, help="Control BED + name")
    pa.add_argument("-d", action="append", nargs=2, metavar=("CSV","NAME"),
                    required=True, help="DE CSV + name (locus_tag, log2FC)")
    pa.add_argument("--operon", required=False, help="2-col TSV mapping operon_id to locus_tag (used only if BED tags not in GFF)")
    pa.add_argument("--gff",    required=True, help="Input GFF file for gene locus_tag strands")
    pa.add_argument("-o","--outdir", required=True, help="Output directory")
    pa.add_argument("-p","--prefix", required=True, help="Filename prefix")
    args = pa.parse_args()

    # Load inputs
    gff_strand = read_gff(args.gff)
    de_map     = {n: read_de(p)  for p,n in args.d}
    T_bed      = {n: read_bed(p) for p,n in args.i}
    C_bed      = {n: read_bed(p) for p,n in args.c}
    trts       = sorted(de_map)

    # Decide mode
    mode = detect_mode(T_bed, C_bed, gff_strand)
    print(f"[INFO] Analysis mode: {mode}")

    if mode == "gene":
        df_m, df_o = build_gene_tables(gff_strand, T_bed, C_bed, de_map, trts)
    else:
        if not args.operon or not os.path.isfile(args.operon):
            raise SystemExit("Operon mode detected but --operon map not provided or file missing.")
        op2genes = read_operon_map(args.operon)
        df_m, df_o = build_operon_tables(T_bed, C_bed, de_map, trts, op2genes, gff_strand)

    # Run analyses
    os.makedirs(args.outdir, exist_ok=True)
    scatter_percentiles_with_ci(df_m, trts, f"{args.prefix}_match", args.outdir)
    scatter_percentiles_with_ci(df_o, trts, f"{args.prefix}_opp",   args.outdir)

if __name__ == "__main__":
    main()
