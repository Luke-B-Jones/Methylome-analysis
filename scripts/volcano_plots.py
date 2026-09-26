#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate volcano plots from multiple DESeq2 output CSV files.'
    )
    parser.add_argument(
        '-i', '--inputs', nargs='+', required=True, metavar=('FILE', 'PREFIX'),
        help='Pairs of <CSV file path> <prefix name for output image>. Can be repeated.'
    )
    parser.add_argument(
        '-o', '--output_dir', required=True,
        help='Output directory to store all plots.'
    )
    parser.add_argument(
        '--padj', type=float, default=0.05,
        help='Adjusted p-value threshold for significance (default: 0.05).'
    )
    parser.add_argument(
        '--log2fc', type=float, default=0.5,
        help='Log2 fold change threshold for significance.'
    )
    parser.add_argument(
        '--x', type=float, help='Set manual x-axis range as [-x, x]'
    )
    parser.add_argument(
        '--y', type=float, help='Set manual y-axis max (y-axis will be [0, y])'
    )
    return parser.parse_args()

def read_data(input_file):
    try:
        df = pd.read_csv(input_file, index_col=0)
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'gene'}, inplace=True)
        required_cols = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
        missing = list(set(required_cols) - set(df.columns))
        if missing:
            print(f"Error: Missing columns in input file {input_file}: {missing}")
            sys.exit(1)
        return df
    except Exception as e:
        print(f"Error reading the input file {input_file}: {e}")
        sys.exit(1)


def prepare_data(df, padj_threshold=0.05, log2fc_threshold=0.5):
    min_padj = df.loc[df['padj'] > 0, 'padj'].min()
    if pd.notnull(min_padj):
        df['padj'] = df['padj'].replace(0, min_padj)
        df.loc[df['padj'] < 0, 'padj'] = min_padj
    else:
        df['padj'] = df['padj'].replace(0, 1e-300)
        df.loc[df['padj'] < 0, 'padj'] = 1e-300

    df['-log10(padj)'] = -np.log10(df['padj'])

    conditions = [
        (df['padj'] < padj_threshold) & (df['log2FoldChange'] > log2fc_threshold),
        (df['padj'] < padj_threshold) & (df['log2FoldChange'] < -log2fc_threshold)
    ]
    choices = ['Upregulated', 'Downregulated']
    df['Significance'] = np.select(conditions, choices, default='Not Significant')

    return df

def compute_global_axis_limits(all_dataframes):
    all_x = pd.concat([df['log2FoldChange'] for df in all_dataframes])
    all_y = pd.concat([df['-log10(padj)'] for df in all_dataframes])
    x_lim = (all_x.min() - 0.5, all_x.max() + 0.5)
    y_lim = (0, all_y.max() + 1)
    return x_lim, y_lim

def plot_volcano(df, output_path, title, padj_threshold, log2fc_threshold, xlim, ylim):
    # Print per-treatment counts of up/down regulated genes
    up_count = (df['Significance'] == 'Upregulated').sum()
    down_count = (df['Significance'] == 'Downregulated').sum()
    print(f"{title}/UP={up_count}/DOWN={down_count}")

    # Count clipped points
    clipped_x = df[(df['log2FoldChange'] < xlim[0]) | (df['log2FoldChange'] > xlim[1])].shape[0]
    clipped_y = df[df['-log10(padj)'] > ylim[1]].shape[0]
    total_clipped = clipped_x + clipped_y

    # Print improved output summary
    if total_clipped > 0:
        print(f"{output_path} - {total_clipped} excluded (X: {clipped_x}, Y: {clipped_y})")
    else:
        print(f"{output_path} - 0 excluded")

    # Identify and print outlier genes
    outlier_mask = (
        (df['log2FoldChange'] < xlim[0]) |
        (df['log2FoldChange'] > xlim[1]) |
        (df['-log10(padj)'] > ylim[1])
    )
    outliers = df[outlier_mask]

    if not outliers.empty:
        print(f"\nOutlier genes in '{title}':")
        for gene in outliers['gene']:
            print(f"  {gene}")

    plt.figure(figsize=(10, 8))

    colors = {'Upregulated': '#0072B2', 'Downregulated': '#db6f1d', 'Not Significant': 'gray'}
    alphas = {'Upregulated': 0.6, 'Downregulated': 0.6, 'Not Significant': 0.25}

    for category in colors:
        subset = df[df['Significance'] == category]
        plt.scatter(
            subset['log2FoldChange'], subset['-log10(padj)'],
            c=colors[category],
            alpha=alphas[category],
            edgecolors='none',
            s=20
        )

    plt.axhline(-np.log10(padj_threshold), color='grey', linestyle='dashed', linewidth=0.75)
    plt.axvline(log2fc_threshold, color='grey', linestyle='dashed', linewidth=0.75)
    plt.axvline(-log2fc_threshold, color='grey', linestyle='dashed', linewidth=0.75)

    plt.xlabel("Log\u2082 Fold Change", fontsize=28)  # \u2082 = subscript 2
    plt.ylabel("-Log\u2081\u2080(p-value)", fontsize=28)  # \u2081\u2080 = subscript 10


    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)

    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.title(title, fontsize=18)
    plt.grid(False)
    plt.tight_layout()

    try:
        plt.savefig(output_path, dpi=300)
    except Exception as e:
        print(f"Error saving plot to {output_path}: {e}")
        sys.exit(1)
    finally:
        plt.close()


def main():
    args = parse_arguments()

    if len(args.inputs) % 2 != 0:
        print("Error: The -i / --inputs argument must contain pairs of <file> <prefix>.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    input_pairs = list(zip(args.inputs[::2], args.inputs[1::2]))

    all_prepared = []
    for filepath, _ in input_pairs:
        if not os.path.isfile(filepath):
            print(f"Error: Input file '{filepath}' does not exist.")
            sys.exit(1)
        df = read_data(filepath)
        prepared_df = prepare_data(df, args.padj, args.log2fc)
        all_prepared.append(prepared_df)

    # Axis limit handling
    if args.x is not None:
        xlim = (-args.x, args.x)
    else:
        xlim, _ = compute_global_axis_limits(all_prepared)

    if args.y is not None:
        ylim = (0, args.y)
    else:
        _, ylim = compute_global_axis_limits(all_prepared)

    for (filepath, prefix), df in zip(input_pairs, all_prepared):
        out_file = os.path.join(args.output_dir, f"{prefix}.png")
        plot_volcano(df, out_file, title=prefix, padj_threshold=args.padj,
                     log2fc_threshold=args.log2fc, xlim=xlim, ylim=ylim)

if __name__ == "__main__":
    main()
