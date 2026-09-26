#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Hypergeometric enrichment of significant bases in feature sets "
            "(promoters/operons/genes/terminators) against a background universe."
        )
    )
    p.add_argument('--full', nargs=2, action='append', metavar=('FILE','SAMPLE'),
                   required=True, help="BED-like of ALL modifiable bases (with header) + sample name")
    p.add_argument('--sig',  nargs=2, action='append', metavar=('FILE','SAMPLE'),
                   required=True, help="BED-like of SIGNIFICANT bases (with header) + sample name")
    p.add_argument('-b','--background', required=True,
                   help="Background regions (BED, with or without header)")
    # Optional feature sets (BED). Provide any subset.
    p.add_argument('--promoters',   help="Promoters BED (upstream regions; optional)")
    p.add_argument('--operons',     help="Operons BED (optional)")
    p.add_argument('--genes',       help="Genes BED (optional)")
    p.add_argument('--terminators', help="Terminators BED (optional)")
    p.add_argument('-o','--outdir', required=True, help="Output directory")
    p.add_argument('-p','--prefix', required=True, help="Output filename prefix")

    # FDR scope: choose the family of hypotheses to control
    p.add_argument('--fdr-scope',
                   choices=['by-sample', 'by-feature', 'by-mark-feature', 'global'],
                   default='by-sample',
                   help=("Family for BH/FDR: "
                         "'by-sample' = features within each sample; "
                         "'by-feature' = samples within each feature; "
                         "'by-mark-feature' = samples within each (prefix, feature); "
                         "'global' = all rows together. "
                         "With m=1, BH returns p_adj=p."))
    return p.parse_args()

def load_bed_points(path: str) -> pd.DataFrame:
    """
    Load a BED/bed-like table WITH header where the first three columns are chrom,start,end.
    Coerce to numeric start/end and return DataFrame with columns ['chrom','start','end'].
    Filters invalid rows (nan start/end or end <= start). Resets index to contiguous [0..n-1].
    """
    df = pd.read_csv(path, sep='\t', header=0, dtype=str)
    cols = list(df.columns)
    if len(cols) < 3:
        raise ValueError(f"{path}: expected at least 3 columns (chrom,start,end)")
    out = pd.DataFrame({
        'chrom': df.iloc[:,0].astype(str),
        'start': pd.to_numeric(df.iloc[:,1], errors='coerce'),
        'end':   pd.to_numeric(df.iloc[:,2], errors='coerce'),
    }).dropna(subset=['start','end'])
    out['start'] = out['start'].astype(int)
    out['end']   = out['end'].astype(int)
    out = out[out['end'] > out['start']].reset_index(drop=True)
    return out

def load_bed_regions(path: str):
    """
    Load a BED file that may or may not have a header.
    Returns list of (chrom,start,end) with start/end as ints; invalid rows are dropped.
    """
    def to_list(df: pd.DataFrame):
        df = df.dropna(subset=[df.columns[1], df.columns[2]])
        chrom = df.iloc[:,0].astype(str)
        start = pd.to_numeric(df.iloc[:,1], errors='coerce')
        end   = pd.to_numeric(df.iloc[:,2], errors='coerce')
        mask = start.notna() & end.notna() & (end.values > start.values)
        return list(zip(chrom[mask].tolist(), start[mask].astype(int).tolist(), end[mask].astype(int).tolist()))

    try:
        df = pd.read_csv(path, sep='\t', header=0, dtype=str)
        regs = to_list(df)
        if len(regs) == 0:
            df = pd.read_csv(path, sep='\t', header=None, dtype=str)
            regs = to_list(df)
    except Exception:
        df = pd.read_csv(path, sep='\t', header=None, dtype=str)
        regs = to_list(df)
    return regs

def build_overlap_mask(df: pd.DataFrame, regs) -> np.ndarray:
    """
    Vectorized interval overlap using positional indices only:
      an interval [s,e) overlaps [rs,re) if (s < re) and (e > rs).
    df: DataFrame with columns chrom,start,end
    regs: list[(chrom,start,end)]
    Returns boolean numpy array mask (len == len(df)).
    """
    n = len(df)
    mask = np.zeros(n, dtype=bool)
    if n == 0 or not regs:
        return mask

    chrom_arr = df['chrom'].to_numpy()
    start_arr = df['start'].to_numpy()
    end_arr   = df['end'].to_numpy()

    chrom_to_pos = {}
    for c in pd.unique(chrom_arr):
        chrom_to_pos[c] = np.where(chrom_arr == c)[0]

    for chrom, rs, re in regs:
        pos_idx = chrom_to_pos.get(chrom)
        if pos_idx is None or pos_idx.size == 0:
            continue
        s = start_arr[pos_idx]
        e = end_arr[pos_idx]
        m = (s < re) & (e > rs)
        if m.any():
            mask[pos_idx[m]] = True
    return mask

def enrich_for_feature(full_df: pd.DataFrame, sig_df: pd.DataFrame, bg_regs, fg_regs):
    """
    Compute (N, K, n, k, p) for one feature set.
    N: #full in background
    K: #full in (feature ∩ background)
    n: #sig  in background
    k: #sig  in (feature ∩ background)
    """
    full_in_bg = build_overlap_mask(full_df, bg_regs)
    sig_in_bg  = build_overlap_mask(sig_df,  bg_regs)

    N = int(full_in_bg.sum())
    n = int(sig_in_bg.sum())
    if N == 0 or n == 0:
        return N, 0, n, 0, 1.0

    fg_full_bg = full_in_bg & build_overlap_mask(full_df, fg_regs)
    fg_sig_bg  = sig_in_bg  & build_overlap_mask(sig_df,  fg_regs)

    K = int(fg_full_bg.sum())
    k = int(fg_sig_bg.sum())
    if K == 0 or k == 0:
        return N, K, n, k, 1.0

    p = float(hypergeom.sf(k - 1, N, K, n))
    return N, K, n, k, p

def _pairs_to_map(pairs, flag_name):
    """
    Build a mapping {sample_name: file_path} from a list of (file, sample) pairs.
    Errors on duplicate sample names to avoid silent overwrites.
    """
    mapping = {}
    for f, s in pairs:
        if s in mapping:
            raise SystemExit(f"Duplicate sample name {s!r} in {flag_name}. Each sample must be unique.")
        mapping[s] = f
    return mapping

def _assign_fdr(df: pd.DataFrame, indexer):
    """Assign BH/FDR to df['p_adj'] for the subset defined by positional/label indexer."""
    pv = df.loc[indexer, 'p_value'].to_numpy(dtype=float)
    if pv.size == 0:
        return
    adj = multipletests(pv, method='fdr_bh')[1]
    df.loc[indexer, 'p_adj'] = adj

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Collect features based on provided files
    features = []
    if args.promoters:
        features.append(("promoters", load_bed_regions(args.promoters)))
    if args.operons:
        features.append(("operons", load_bed_regions(args.operons)))
    if args.genes:
        features.append(("genes", load_bed_regions(args.genes)))
    if args.terminators:
        features.append(("terminators", load_bed_regions(args.terminators)))

    if len(features) == 0:
        raise SystemExit("No features provided. Supply at least one of --promoters/--operons/--genes/--terminators.")

    bg_regs = load_bed_regions(args.background)

    # Map sample → file path with duplicate-sample protection
    full_map = _pairs_to_map(args.full, "--full")
    sig_map  = _pairs_to_map(args.sig,  "--sig")

    if set(full_map) != set(sig_map):
        missing = (set(full_map) ^ set(sig_map))
        raise SystemExit(f"Sample mismatch between --full and --sig: {sorted(missing)}")

    rows = []

    # Per-sample calculations
    for sample in sorted(full_map.keys()):
        full_df = load_bed_points(full_map[sample])  # header present
        sig_df  = load_bed_points(sig_map[sample])   # header present

        for label, fg_regs in features:
            N, K, n, k, p = enrich_for_feature(full_df, sig_df, bg_regs, fg_regs)
            # Only keep informative rows (requires at least one sig in bg and in feature∩bg)
            if n == 0 or k == 0:
                continue
            rows.append({
                'Sample':  sample,
                'Feature': label,
                'N': N, 'K': K, 'n': n, 'k': k,
                'p_value': p
            })

    # Combined "ALL" across samples; ensure positional indices are contiguous
    all_full = (
        pd.concat([load_bed_points(full_map[s]) for s in full_map], ignore_index=True)
          .drop_duplicates()
          .reset_index(drop=True)
    )
    all_sig = (
        pd.concat([load_bed_points(sig_map[s]) for s in sig_map], ignore_index=True)
          .drop_duplicates()
          .reset_index(drop=True)
    )

    for label, fg_regs in features:
        N, K, n, k, p = enrich_for_feature(all_full, all_sig, bg_regs, fg_regs)
        if n == 0 or k == 0:
            continue
        rows.append({
            'Sample':  'ALL',
            'Feature': label,
            'N': N, 'K': K, 'n': n, 'k': k,
            'p_value': p
        })

    # Output table (+ configurable FDR)
    df = pd.DataFrame(rows)
    if not df.empty:
        df['p_adj'] = np.nan

        scope = args.fdr_scope
        if scope == 'by-sample':
            for _, idx in df.groupby('Sample').groups.items():
                _assign_fdr(df, idx)

        elif scope == 'by-feature':
            for _, idx in df.groupby('Feature').groups.items():
                _assign_fdr(df, idx)

        elif scope == 'by-mark-feature':
            # 'mark' inferred from this run's prefix (keeps A vs C universes separate when you mix runs)
            df['_mark'] = args.prefix
            for _, idx in df.groupby(['_mark', 'Feature']).groups.items():
                _assign_fdr(df, idx)
            df.drop(columns=['_mark'], inplace=True)

        elif scope == 'global':
            _assign_fdr(df, df.index)
    else:
        df = pd.DataFrame(columns=['Sample','Feature','N','K','n','k','p_value','p_adj'])

    out_path = os.path.join(args.outdir, f"{args.prefix}_hypergeom.tsv")
    df.to_csv(out_path, sep='\t', index=False)
    # Console echo for convenience
    print(df.to_string(index=False))
    print("Written:", out_path)

if __name__ == "__main__":
    main()
