#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    p = argparse.ArgumentParser(
        description="Map DESeq2 log2FC to KOs and aggregate by KO (mean)."
    )
    p.add_argument('--deseq2', required=True,
                   help="DESeq2 CSV with first column = locus tag, and a 'log2FoldChange' column")
    p.add_argument('--ko', required=True,
                   help="KO mapping TSV, two columns: locus_tag[TAB]KO (no header).")
    p.add_argument('-o', '--output', required=True,
                   help="Output TSV: KO<TAB>log2FoldChange (mean across genes).")
    args = p.parse_args()

    # 1) Read DESeq2; first column is locus tag
    try:
        df = pd.read_csv(args.deseq2, index_col=0)
    except Exception as e:
        sys.stderr.write(f"Error reading DESeq2 file: {e}\n")
        sys.exit(1)
    if 'log2FoldChange' not in df.columns:
        sys.stderr.write("Error: 'log2FoldChange' column not found in DESeq2 file.\n")
        sys.exit(1)

    # 2) Read KO mapping
    try:
        ko_df = pd.read_csv(
            args.ko,
            sep='\t',
            header=None,
            names=['locus_tag', 'KO'],
            dtype={'locus_tag': str, 'KO': str},
            keep_default_na=False
        )
    except Exception as e:
        sys.stderr.write(f"Error reading KO mapping: {e}\n")
        sys.exit(1)

    # Drop any rows with empty KO
    ko_df = ko_df[ko_df['KO'].str.strip() != ""]

    # 3) Merge on locus_tag
    #    Reset index of DESeq2 to bring locus_tag back as a column
    df2 = df.reset_index().rename(columns={'index':'locus_tag'})
    merged = pd.merge(
        ko_df, 
        df2[['locus_tag','log2FoldChange']],
        on='locus_tag',
        how='inner'
    )

    if merged.empty:
        sys.stderr.write("Warning: No overlapping locus_tags found between DESeq2 and KO mapping.\n")

    # 4) Aggregate: mean log2FoldChange per KO, skipping NaNs
    agg = ( merged
            .groupby('KO', as_index=False)
            ['log2FoldChange']
            .mean()
            .dropna(subset=['log2FoldChange'])
          )

    # 5) Write output
    agg.to_csv(
        args.output,
        sep='\t',
        index=False,
        float_format="%.6f"
    )
    print(f"Wrote {len(agg)} KOs to {args.output}")

if __name__ == '__main__':
    main()
