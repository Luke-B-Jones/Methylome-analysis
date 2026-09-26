#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

# Transparency levels for plotting layers
ALL_ALPHA = 0.5
KO_ALPHA = 0.6

# Color, transparency, and size for 'x' markers
X_COLOR = 'black'
X_ALPHA = 0.5
X_SIZE = 15  # size of the 'x' markers

# Define colors for each treatment: (up_color, down_color)
TREATMENT_COLORS = {
    'H2O2': ('#1b9e77', '#8dcfbb'),
    'AMPI': ('#d95f02', '#ec8780'),
    'VANC': ('#7570b3', '#bab8d9'),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a ridge plot of Log2 Fold Change distributions per treatment. "
            "First layer: all genes in grey; second layer: only KO‐assigned genes in color. "
            "Additionally, mark significant KO genes (|Log2FC| > 0.5 and adjusted p‐value < 0.05) "
            "with a small 'x', and annotate counts of KO genes on each side of Log2FC = 0."
        )
    )
    parser.add_argument(
        '-i', '--index',
        required=True,
        help='Path to the index file (TSV) mapping locus tags to KO assignments.'
    )
    parser.add_argument(
        '--de',
        action='append',
        nargs=2,
        metavar=('DE_FILE', 'TREATMENT_NAME'),
        required=True,
        help=(
            "Specify a differential‐expression CSV file and its treatment name. "
            "Usage: --de <path/to/DE.csv> <TREATMENT_NAME>. Can be given multiple times."
        )
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to the output image (e.g., output.png).'
    )
    return parser.parse_args()


def load_index(index_path):
    """
    Load the index file (TSV with two columns: gene identifier and KO assignment).
    Extract the locus tag (last component after '|') and keep only those with a KO.
    Returns a set of locus tags that have a KO assignment.
    """
    idx_df = pd.read_csv(index_path, sep='\t', header=None, names=['gene', 'KO'], dtype=str)
    idx_df['locus'] = idx_df['gene'].apply(lambda x: x.split('|')[-1] if isinstance(x, str) else '')
    has_ko = idx_df['KO'].notna() & (idx_df['KO'].str.strip() != '')
    loci_with_ko = set(idx_df.loc[has_ko, 'locus'])
    return loci_with_ko


def compute_kde(values, x_grid):
    """
    Given a 1D array of values, compute the Gaussian KDE over x_grid.
    Returns a density array (same shape as x_grid), normalized so that its
    maximum is 1.
    """
    kde = gaussian_kde(values)
    density = kde(x_grid)
    if density.max() > 0:
        density = density / density.max()
    return density


def main():
    args = parse_args()

    loci_with_ko = load_index(args.index)
    if not loci_with_ko:
        sys.exit("ERROR: No loci with KO assignments found in the index file.")

    treatment_data = {}
    treatment_order = []

    # Thresholds for significance
    LOG2FC_THRESHOLD = 0.5
    ADJ_PVALUE_THRESHOLD = 0.05

    for de_file, treatment_name in args.de:
        df = pd.read_csv(de_file, index_col=0)
        df.index = df.index.map(lambda x: x.strip('"').strip("'") if isinstance(x, str) else x)
        if 'log2FoldChange' not in df.columns or 'padj' not in df.columns:
            sys.exit(f"ERROR: 'log2FoldChange' and 'padj' columns must be present in {de_file}")

        all_vals = df['log2FoldChange'].dropna().values.astype(float)
        ko_mask = df.index.isin(loci_with_ko)
        ko_vals = df.loc[ko_mask, 'log2FoldChange'].dropna().values.astype(float)

        # Identify significant KO genes based on adjusted p‐value
        sig_mask = ko_mask & (df['log2FoldChange'].abs() > LOG2FC_THRESHOLD) & (df['padj'] < ADJ_PVALUE_THRESHOLD)
        sig_vals = df.loc[sig_mask, 'log2FoldChange'].dropna().values.astype(float)

        if all_vals.size == 0:
            print(f"Warning: No valid log2FoldChange values in DE file '{de_file}'. Skipping.")
            continue

        treatment_data[treatment_name] = {
            'all': all_vals,
            'ko': ko_vals,
            'sig': sig_vals
        }
        treatment_order.append(treatment_name)

    if not treatment_data:
        sys.exit("ERROR: No valid DE data loaded for any treatment.")

    # Define x‐axis range capped at ±4
    x_min, x_max = -4, 4
    x_grid = np.linspace(x_min, x_max, 1000)  # High resolution for KDE

    num_treatments = len(treatment_order)
    scale = 1.0
    buffer = 0.6  # Slightly increased buffer for spacing
    y_positions = np.arange(num_treatments) * (scale + buffer)

    fig_height = num_treatments * (scale + buffer) + 1.0
    fig, ax = plt.subplots(figsize=(8, fig_height))

    epsilon = (x_grid[1] - x_grid[0]) / 2

    # Fixed horizontal positions for annotation text (closer to center)
    left_x_pos = -1.0
    right_x_pos = 1.0

    # Fixed vertical offset from the ridge peak for annotation text
    y_text_offset = scale + buffer * 0.1

    for idx, treatment in enumerate(treatment_order):
        data = treatment_data[treatment]
        all_vals = data['all']
        ko_vals = data['ko']
        sig_vals = data['sig']

        # Clip values to x‐axis range for KDE computation
        all_vals_clipped = all_vals[(all_vals >= x_min) & (all_vals <= x_max)]
        density_all = compute_kde(all_vals_clipped, x_grid)
        d_all_scaled = density_all * scale

        y_base = y_positions[idx]

        # Plot grey layer for all genes
        ax.fill_between(
            x_grid,
            y_base,
            y_base + d_all_scaled,
            color='#aaaaaa',
            alpha=ALL_ALPHA,
            linewidth=0
        )

        if ko_vals.size > 0:
            ko_vals_clipped = ko_vals[(ko_vals >= x_min) & (ko_vals <= x_max)]
            density_ko = compute_kde(ko_vals_clipped, x_grid)
            d_ko_scaled = density_ko * scale

            color_up, color_down = TREATMENT_COLORS.get(
                treatment,
                ('#333333', '#cccccc')
            )

            # Fill negative side
            ax.fill_between(
                x_grid,
                y_base,
                y_base + d_ko_scaled,
                where=(x_grid <= epsilon),
                color=color_down,
                alpha=KO_ALPHA,
                linewidth=0
            )
            # Fill positive side
            ax.fill_between(
                x_grid,
                y_base,
                y_base + d_ko_scaled,
                where=(x_grid >= -epsilon),
                color=color_up,
                alpha=KO_ALPHA,
                linewidth=0
            )
            # Outline ridge
            ax.plot(x_grid, y_base + d_ko_scaled, color='black', linewidth=0.6)

        # Plot significant KO gene points as tiny 'x' markers
        if sig_vals.size > 0:
            sig_vals_clipped = sig_vals[(sig_vals >= x_min) & (sig_vals <= x_max)]
            ax.scatter(
                sig_vals_clipped,
                [y_base] * len(sig_vals_clipped),
                marker='x',
                color=X_COLOR,
                alpha=X_ALPHA,
                s=X_SIZE,
                linewidths=0.5
            )

        # Compute counts for annotation
        pos_all_count = np.sum(all_vals_clipped > 0)
        neg_all_count = np.sum(all_vals_clipped < 0)

        pos_ko_count = np.sum(ko_vals_clipped > 0) if ko_vals.size > 0 else 0
        neg_ko_count = np.sum(ko_vals_clipped < 0) if ko_vals.size > 0 else 0

        pos_percent = (pos_ko_count / pos_all_count * 100) if pos_all_count > 0 else 0
        neg_percent = (neg_ko_count / neg_all_count * 100) if neg_all_count > 0 else 0

        # Prepare annotation text (no backslash before %)
        left_text = f"$\\mathit{{n}}$ = {neg_ko_count} ({neg_percent:.1f}%)"
        right_text = f"$\\mathit{{n}}$ = {pos_ko_count} ({pos_percent:.1f}%)"

        # Uniform vertical position for annotation text above each ridge
        y_text = y_base + y_text_offset

        # Left side annotation
        ax.text(
            left_x_pos,
            y_text,
            left_text,
            ha='right',
            va='center',
            fontsize=10,
            fontfamily='sans-serif',
            color='#555555'
        )
        # Right side annotation
        ax.text(
            right_x_pos,
            y_text,
            right_text,
            ha='left',
            va='center',
            fontsize=10,
            fontfamily='sans-serif',
            color='#555555'
        )

    # Add a vertical dashed line at x = 0
    ax.axvline(
        0,
        color='#555555',
        linestyle='--',
        linewidth=0.6
    )

    # Set x‐axis label
    ax.set_xlabel(r'Log$_2$ Fold Change', fontsize=14, fontfamily='sans-serif')

    # Set y‐axis ticks and labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        treatment_order,
        fontsize=12,
        fontfamily='sans-serif'
    )
    ax.tick_params(axis='y', which='both', length=4, width=1)

    # Set axis limits
    ax.set_ylim(-0.5, y_positions[-1] + scale + buffer + 0.2)
    ax.set_xlim(x_min, x_max)

    # Tighten layout and save with higher DPI
    plt.tight_layout()
    plt.savefig(args.output, dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
