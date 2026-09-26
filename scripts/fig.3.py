#!/usr/bin/env python3
import argparse
import concurrent.futures
import hashlib
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
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

# Absolute diff -> alpha mapping (unchanged)
ALPHA_MAP = {(0.2, 0.4): 0.2, (0.4, 0.6): 0.4, (0.6, 0.8): 0.8}
ALPHA_MAX = 1.0

# ----------------------------
# Scaling + DPI (unchanged)
# ----------------------------
FONT_SCALE = 1.2
DPI = 1600

# ----------------------------
# Supported mods
# ----------------------------
VALID_MODS = ('4mC', '5mC', '6mA')


def normalize_mod(mod):
    """Normalize user-supplied mod labels to one of: 4mC, 5mC, 6mA."""
    lookup = {
        '4mc': '4mC',
        '5mc': '5mC',
        '6ma': '6mA',
    }
    key = mod.strip().lower()
    if key not in lookup:
        raise ValueError(
            f"Unsupported mod '{mod}'. Valid mods are: {', '.join(VALID_MODS)}"
        )
    return lookup[key]


def normalize_color_pair(color_pair):
    """
    Parse and validate a user-supplied colour pair like:
    '#1b9e77,#8dcfbb'
    Returns a normalized 2-tuple of hex colours.
    """
    parts = [p.strip().strip('"').strip("'") for p in color_pair.split(',')]
    if len(parts) != 2 or not all(parts):
        raise ValueError(
            f"Invalid colour pair '{color_pair}'. Expected format: '#1b9e77,#8dcfbb'"
        )

    try:
        c1 = mcolors.to_hex(mcolors.to_rgba(parts[0]), keep_alpha=False)
        c2 = mcolors.to_hex(mcolors.to_rgba(parts[1]), keep_alpha=False)
    except ValueError as e:
        raise ValueError(
            f"Invalid colour pair '{color_pair}'. Expected format: '#1b9e77,#8dcfbb'"
        ) from e

    return c1, c2


def build_sample_order(input_specs):
    """Preserve sample order by first appearance on the command line."""
    seen = set()
    order = []
    for _, sample, _, _ in input_specs:
        if sample not in seen:
            seen.add(sample)
            order.append(sample)
    return order


def build_sample_colors(input_specs):
    """
    Build per-sample colour mapping from CLI input.
    Enforces that repeated uses of the same sample keep the same colour pair.
    """
    sample_colors = {}
    for _, sample, _, color_pair in input_specs:
        if sample in sample_colors:
            if sample_colors[sample] != color_pair:
                raise ValueError(
                    f"Sample '{sample}' was given conflicting colour pairs: "
                    f"{sample_colors[sample]} vs {color_pair}"
                )
        else:
            sample_colors[sample] = color_pair
    return sample_colors


def load_and_prepare(path, sample, mod):
    df = pd.read_csv(
        path,
        sep='\t',
        header=None,
        skiprows=1,
        names=COL_NAMES,
        dtype=DTYPE_DICT,
        low_memory=False
    )
    df['control'] = pd.to_numeric(df['control_pct_modified'], errors='coerce') * 100
    df['treatment'] = pd.to_numeric(df['treatment_pct_modified'], errors='coerce') * 100
    df['diff'] = pd.to_numeric(df['differential_pct_modification'], errors='coerce')
    df['sample'] = sample
    df['mod'] = mod
    return df[['sample', 'mod', 'control', 'treatment', 'diff']]


def compute_alpha(d):
    if pd.isna(d) or abs(d) < 0.2:
        return None
    for (lo, hi), a in ALPHA_MAP.items():
        if lo <= abs(d) < hi:
            return a
    return ALPHA_MAX


def save_legend_image(handles, labels, outfile, ncol=1, fontsize=5.5):
    """Save legend only as its own image (key). Returns True if written."""
    if not handles:
        return False

    fig = plt.figure(figsize=(4, 0.6 + 0.3 * len(labels)))
    fig.legend(
        handles,
        labels,
        loc='center',
        frameon=False,
        ncol=ncol,
        fontsize=fontsize * FONT_SCALE
    )
    fig.savefig(outfile, dpi=DPI, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    return True


def plot_mod(df, mod, base, ext, out_folder, sample_order, sample_colors):
    fig, ax = plt.subplots(figsize=(3, 4))
    subset = df[df['mod'] == mod].copy()
    if subset.empty:
        plt.close(fig)
        return

    subset = subset.assign(x0=X_CTRL, x1=X_TRMT)

    handles, labels = [], []

    for sample in sample_order:
        gdf = subset[subset['sample'] == sample].copy()
        if gdf.empty:
            continue

        gdf['alpha'] = gdf['diff'].apply(compute_alpha)
        gdf.dropna(subset=['alpha'], inplace=True)
        if gdf.empty:
            continue

        line_color = sample_colors[sample][0]

        for a, sub in gdf.groupby('alpha'):
            segs = [
                ((x0, y0), (x1, y1))
                for x0, y0, x1, y1 in zip(sub['x0'], sub['control'], sub['x1'], sub['treatment'])
            ]
            lc = LineCollection(
                segs,
                linewidths=LINE_WIDTH,
                colors=line_color,
                alpha=a
            )
            ax.add_collection(lc)

        hypo = (gdf['diff'] < 0).sum()
        hyper = (gdf['diff'] > 0).sum()
        if hypo + hyper:
            handles.append(Line2D([], [], color=line_color, linewidth=2))
            labels.append(f"{sample}: ↑ {hyper}, ↓ {hypo} (ΔPL)")

    ax.set_xticks([X_CTRL, X_TRMT])
    ax.set_xticklabels(['Pre-treatment', 'Post-treatment'])
    ax.set_xlim(X_CTRL, X_TRMT)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Population-level prevalence (%)', fontsize=8 * FONT_SCALE)

    ax.tick_params(axis='x', labelsize=6.5 * FONT_SCALE)
    ax.tick_params(axis='y', labelsize=6.5 * FONT_SCALE)

    if handles:
        ax.legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.075),
            frameon=False,
            fontsize=5.5 * FONT_SCALE,
            ncol=1
        )
        fig.subplots_adjust(bottom=0.3)

    fig.tight_layout()
    outfile = os.path.join(out_folder, f"{base}_{mod}{ext}")
    fig.savefig(outfile, dpi=DPI)
    plt.close(fig)
    print(f"Saved plot: {outfile}")

    legend_out = os.path.join(out_folder, f"{base}_{mod}_legend{ext}")
    if save_legend_image(handles, labels, legend_out, ncol=1):
        print(f"Saved legend: {legend_out}")


def make_cache_signature(input_specs):
    """
    Build a cache signature from file paths, sizes, mtimes, samples, mods and colours.
    This keeps --resume safe when flexible -i inputs change.
    """
    parts = []
    for path, sample, mod, color_pair in input_specs:
        st = os.stat(path)
        parts.append(
            '\t'.join([
                os.path.abspath(path),
                sample,
                mod,
                color_pair[0],
                color_pair[1],
                str(st.st_size),
                str(st.st_mtime_ns),
            ])
        )
    return hashlib.sha256('\n'.join(parts).encode('utf-8')).hexdigest()


def parse_inputs(raw_inputs):
    input_specs = []
    for path, sample, mod, color_pair in raw_inputs:
        norm_mod = normalize_mod(mod)
        norm_colors = normalize_color_pair(color_pair)
        input_specs.append((path, sample, norm_mod, norm_colors))
    return input_specs


def main():
    parser = argparse.ArgumentParser(
        description="Paired line plots with transparency for population-level prevalence."
    )
    parser.add_argument('-png', required=True, help='Base output filename, e.g. summary.png')
    parser.add_argument('-o', required=True, help='Output directory')
    parser.add_argument('--resume', action='store_true', help='Reuse cached combined dataframe if inputs match')
    parser.add_argument(
        '-i',
        action='append',
        nargs=4,
        metavar=('FILE', 'SAMPLE', 'MOD', 'COLORS'),
        required=True,
        help=(
            "Add one input as: -i <file.bed> <sample> <mod> <color1,color2>. "
            "Repeat as needed. MOD must be one of: 4mC, 5mC, 6mA. "
            "Example colours: '#1b9e77,#8dcfbb'"
        )
    )

    args = parser.parse_args()

    try:
        input_specs = parse_inputs(args.i)
        sample_colors = build_sample_colors(input_specs)
    except ValueError as e:
        parser.error(str(e))

    for path, _, _, _ in input_specs:
        if not os.path.exists(path):
            parser.error(f"Input file not found: {path}")

    out = args.o
    os.makedirs(out, exist_ok=True)

    sample_order = build_sample_order(input_specs)

    cache = os.path.join(out, 'cached.pkl')
    cache_meta = os.path.join(out, 'cached.meta.json')
    signature = make_cache_signature(input_specs)

    use_cache = False
    if args.resume and os.path.exists(cache) and os.path.exists(cache_meta):
        try:
            with open(cache_meta, 'r', encoding='utf-8') as fh:
                meta = json.load(fh)
            use_cache = (meta.get('signature') == signature)
        except Exception:
            use_cache = False

    if use_cache:
        df_all = pd.read_pickle(cache)
        print(f"Loaded dataset from {cache}")
    else:
        dfs = []
        with concurrent.futures.ThreadPoolExecutor() as ex:
            futures = {
                ex.submit(load_and_prepare, path, sample, mod): (path, sample, mod)
                for path, sample, mod, _ in input_specs
            }
            for fut, (path, sample, mod) in futures.items():
                try:
                    dfs.append(fut.result())
                except Exception as e:
                    print(f"Error {sample} / {mod} / {path}: {e}")

        if not dfs:
            raise RuntimeError("No valid input datasets could be loaded.")

        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_pickle(cache)
        with open(cache_meta, 'w', encoding='utf-8') as fh:
            json.dump({'signature': signature}, fh)
        print(f"Saved combined dataset to {cache}")

    base, ext = os.path.splitext(args.png)
    for mod in VALID_MODS:
        plot_mod(df_all, mod, base, ext, out, sample_order, sample_colors)


if __name__ == '__main__':
    main()