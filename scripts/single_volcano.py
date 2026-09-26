#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate a volcano plot from DESeq2 output, then print up/down lists.'
    )
    parser.add_argument(
        '-i', '--input', required=True, help='Path to DESeq2 output CSV file.'
    )
    parser.add_argument(
        '-o', '--output', required=True, 
        help='Prefix for the plot: this will be used as the plot title and to save the file as PREFIX.png in the current directory.'
    )
    parser.add_argument(
        '--padj', type=float, default=0.05, 
        help='Adjusted p-value threshold for significance (default: 0.05).'
    )
    parser.add_argument(
        '--log2fc', type=float, required=True, 
        help='Log2 fold change threshold for significance.'
    )
    return parser.parse_args()

def read_data(input_file):
    try:
        df = pd.read_csv(input_file)
        # Ensure necessary columns are present
        required_cols = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
        missing = list(set(required_cols) - set(df.columns))
        if missing:
            print(f"Error: Missing columns in input file: {missing}")
            sys.exit(1)
        return df
    except Exception as e:
        print(f"Error reading the input file: {e}")
        sys.exit(1)

def prepare_data(df, padj_threshold=0.05, log2fc_threshold=0.5):
    # Replace zero or negative padj with the smallest positive value to avoid -log10 issues
    min_padj = df.loc[df['padj'] > 0, 'padj'].min()
    if pd.notnull(min_padj):
        df['padj'] = df['padj'].replace(0, min_padj)
        df.loc[df['padj'] < 0, 'padj'] = min_padj
    else:
        # Fallback if all padj are zero or NaN
        df['padj'] = df['padj'].replace(0, 1e-300)
        df.loc[df['padj'] < 0, 'padj'] = 1e-300

    df['-log10(padj)'] = -np.log10(df['padj'])
    
    # Define conditions
    conditions = [
        (df['padj'] < padj_threshold) & (df['log2FoldChange'] > log2fc_threshold),
        (df['padj'] < padj_threshold) & (df['log2FoldChange'] < -log2fc_threshold)
    ]
    choices = ['Upregulated', 'Downregulated']
    df['Significance'] = np.select(conditions, choices, default='Not Significant')
    
    return df

def plot_volcano(df, prefix, padj_threshold=0.05, log2fc_threshold=0.5):
    plt.figure(figsize=(10,8))
    
    # Define colors and alpha values for different gene groups
    colors = {'Upregulated': 'green', 'Downregulated': 'red', 'Not Significant': 'gray'}
    alphas = {'Upregulated': 0.6, 'Downregulated': 0.6, 'Not Significant': 0.25}
    
    # Plot each group without adding a legend
    for category in colors:
        subset = df[df['Significance'] == category]
        plt.scatter(
            subset['log2FoldChange'], subset['-log10(padj)'],
            c=colors[category],
            alpha=alphas[category],
            edgecolors='none',
            s=20
        )
    
    # Add threshold lines (vertical and horizontal)
    plt.axhline(-np.log10(padj_threshold), color='grey', linestyle='dashed', linewidth=0.75)
    plt.axvline(log2fc_threshold, color='grey', linestyle='dashed', linewidth=0.75)
    plt.axvline(-log2fc_threshold, color='grey', linestyle='dashed', linewidth=0.75)
    
    # Set axis labels (titles remain unchanged)
    plt.xlabel(r'Log$_2$ Fold Change', fontsize=14)
    plt.ylabel(r'-Log$_{10}(p\text{-value})$', fontsize=14)
    
    # Increase tick (axis number) font sizes
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Set plot title as the prefix provided
    plt.title(prefix, fontsize=16)
    
    # Remove gridlines for a cleaner look and adjust layout
    plt.grid(False)
    plt.tight_layout()
    
    # Save the plot as PREFIX.png in the current directory
    output_filename = f"{prefix}.png"
    try:
        plt.savefig(output_filename, dpi=300)
        print(f"Volcano plot saved to {output_filename}")
    except Exception as e:
        print(f"Error saving the plot: {e}")
        sys.exit(1)
    finally:
        plt.close()

def main():
    args = parse_arguments()
    
    # Check if the input file exists
    if not os.path.isfile(args.input):
        print(f"Error: The input file {args.input} does not exist.")
        sys.exit(1)
    
    # Read and prepare data
    df = read_data(args.input)
    df_prepared = prepare_data(df, padj_threshold=args.padj, log2fc_threshold=args.log2fc)
    
    # Plot and save using the provided prefix
    plot_volcano(df_prepared, args.output, padj_threshold=args.padj, log2fc_threshold=args.log2fc)
    
    # After plotting, print separate lists for Upregulated and Downregulated genes.
    # The script assumes the first column of the CSV is the gene ID.
    gene_col = df_prepared.columns[0]

    up_list = df_prepared.loc[df_prepared['Significance'] == 'Upregulated', gene_col]
    down_list = df_prepared.loc[df_prepared['Significance'] == 'Downregulated', gene_col]

    print("\nUP-REGULATED GENES:")
    for gene in up_list:
        print(gene)

    print("\nDOWN-REGULATED GENES:")
    for gene in down_list:
        print(gene)

if __name__ == "__main__":
    main()
