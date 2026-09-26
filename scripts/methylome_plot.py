#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import concurrent.futures
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# ----------------------------
# Column definitions and dtypes
# ----------------------------
COL_NAMES = [
    'chrom', 'start', 'end', 'name', 'score', 'strand',
    'control_counts', 'control_total', 'treatment_counts', 'treatment_total',
    'control_mod_percentages', 'treatment_mod_percentages',
    'control_pct_modified', 'treatment_pct_modified', 'differential_pct_modification'
]

DTYPE_DICT = {i: typ for i, typ in enumerate([
    str, int, int, str, float, str,
    str, str, str, str,
    str, str, float, float, float
])}

# ----------------------------
# Plot parameters (unchanged)
# ----------------------------
LINE_WIDTH = 0.4
X_CTRL, X_TRMT = 1.25, 1.75

# Group colors (unchanged)
TREATMENT_COLORS = {
    'H2O2': ('#1b9e77', '#8dcfbb'),
    'AMPI': ('#d95f02', '#ec8780'),
    'VANC': ('#7570b3', '#bab8d9'),
}


# Alpha thresholds for diff (abs) (unchanged)
ALPHA_MAP = {(0.2,0.4):0.2,(0.4,0.6):0.4,(0.6,0.8):0.8}
ALPHA_MAX = 1.0

# ----------------------------
# NEW: scaling + dpi
# ----------------------------
FONT_SCALE = 1.2          # (2) increase font size by 120%
DPI = 1600                # (3) double the resolution (was 800)

def load_and_prepare(path, label):
    df = pd.read_csv(path, sep='\t', header=None, skiprows=1,
                     names=COL_NAMES, dtype=DTYPE_DICT, low_memory=False)
    df['control'] = pd.to_numeric(df['control_pct_modified'], errors='coerce') * 100
    df['treatment'] = pd.to_numeric(df['treatment_pct_modified'], errors='coerce') * 100
    df['diff'] = pd.to_numeric(df['differential_pct_modification'], errors='coerce')
    df['group'] = label
    return df[['group','control','treatment','diff']]

def compute_alpha(d):
    if pd.isna(d) or abs(d) < 0.2:
        return None
    for (lo, hi), a in ALPHA_MAP.items():
        if lo <= abs(d) < hi:
            return a
    return ALPHA_MAX

# NEW: save a separate legend image per plot
def save_legend_image(handles, labels, outfile, ncol=1, fontsize=5.5):
    """Save legend only as its own image (key)."""
    if not handles:
        return
    fig = plt.figure(figsize=(4, 0.6 + 0.3*len(labels)))
    fig.legend(handles, labels, loc='center', frameon=False,
               ncol=ncol, fontsize=fontsize * FONT_SCALE)
    fig.savefig(outfile, dpi=DPI, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

def plot_condition(df, condition, base, ext, out_folder):
    fig, ax = plt.subplots(figsize=(3,4))
    subset = df[df['group'].str.endswith(f"_{condition}")]
    if subset.empty:
        return

    # assign x coords
    subset = subset.assign(x0=X_CTRL, x1=X_TRMT)

    handles, labels = [], []

    for grp in ['H2O2', 'AMPI', 'VANC']:
        grp_label = f"{grp}_{condition}"
        gdf = subset[subset['group'] == grp_label].copy()
        gdf['alpha'] = gdf['diff'].apply(compute_alpha)
        gdf.dropna(subset=['alpha'], inplace=True)
        if gdf.empty:
            continue

        # draw lines per alpha level
        for a, sub in gdf.groupby('alpha'):
            segs = [((x0, y0), (x1, y1)) for x0, y0, x1, y1 in
                    zip(sub['x0'], sub['control'], sub['x1'], sub['treatment'])]
            lc = LineCollection(segs, linewidths=LINE_WIDTH,
                                colors=TREATMENT_COLORS[grp][0], alpha=a)
            ax.add_collection(lc)

        # legend entries
        hypo = (gdf['diff'] < 0).sum()
        hyper = (gdf['diff'] > 0).sum()
        if hypo + hyper:
            handles.append(Line2D([], [], color=TREATMENT_COLORS[grp][0], linewidth=2))
            labels.append(f"{grp}: ↑ {hyper}, ↓ {hypo} (ΔPL)")

    ax.set_xticks([X_CTRL, X_TRMT])
    ax.set_xticklabels(['Pre-treatment', 'Post-treatment'])
    ax.set_xlim(X_CTRL, X_TRMT)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Population-level prevalence (%)', fontsize=8 * FONT_SCALE)

    # (2) bump tick label sizes by 120%
    ax.tick_params(axis='x', labelsize=6.5 * FONT_SCALE)
    ax.tick_params(axis='y', labelsize=6.5 * FONT_SCALE)

    # Keep the on-plot legend (no content lost), but larger fonts
    if handles:
        ax.legend(handles, labels, loc='upper center',
                  bbox_to_anchor=(0.5, -0.075), frameon=False,
                  fontsize=5.5 * FONT_SCALE, ncol=1)
        fig.subplots_adjust(bottom=0.3)

    fig.tight_layout()
    outfile = os.path.join(out_folder, f"{base}_{condition}{ext}")
    fig.savefig(outfile, dpi=DPI)
    plt.close(fig)
    print(f"Saved plot: {outfile}")

    # (1) save a separate legend image (key) for this plot
    legend_out = os.path.join(out_folder, f"{base}_{condition}_legend{ext}")
    save_legend_image(handles, labels, legend_out, ncol=1)
    print(f"Saved legend: {legend_out}")

def main():
    parser = argparse.ArgumentParser(description="Scatter plot with transparency.")
    parser.add_argument('-png', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('--resume', action='store_true')
    for name in ['H2O2_A','H2O2_C','AMPI_A','AMPI_C','VANC_A','VANC_C']:
        parser.add_argument(f'-{name}', required=True)
    args = parser.parse_args()

    out = args.o
    os.makedirs(out, exist_ok=True)
    tasks = [(g, c, getattr(args, f'{g}_{c}'))
             for g in ['H2O2', 'AMPI', 'VANC'] for c in ['A', 'C']]

    cache = os.path.join(out, 'cached.pkl')
    if args.resume and os.path.exists(cache):
        df_all = pd.read_pickle(cache)
        print(f"Loaded dataset from {cache}")
    else:
        dfs = []
        with concurrent.futures.ThreadPoolExecutor() as ex:
            futures = {ex.submit(load_and_prepare, p, label): label
                       for g, c, p in tasks for label in [f"{g}_{c}"]}
            for fut, label in futures.items():
                try:
                    dfs.append(fut.result())
                except Exception as e:
                    print(f"Error {label}: {e}")
        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_pickle(cache)
        print(f"Saved combined dataset to {cache}")

    base, ext = os.path.splitext(args.png)
    for cond in ['A', 'C']:
        plot_condition(df_all, cond, base, ext, out)


if __name__=='__main__':
    main()
