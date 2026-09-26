#!/usr/bin/env python3
"""
map_mods_to_regions.py

Map per-treatment & per-modification-type base modifications onto genome regions,
annotate with product, integrate differential stats, and build a clean tabular report
where positive diffs are prefixed with '+' — and never discard any region that overlaps a mod.

Updated:
- Enforce consistent per-treatment ordering for 4mC, 5mC, 6mA, and Log2FC columns.
- Output a single combined raw TSV containing all modification types.
- Additionally output filtered TSVs split by modification type:
  * <base>_4mC_filtered.tsv
  * <base>_5mC_filtered.tsv
  * <base>_6mA_filtered.tsv
- Filtering removes rows that:
  (a) lack any Log2FC value,
  (b) are gene; with a "hypothetical / ribosomal RNA / tRNA" function, or
  (c) are upstream; where all function entries are "hypothetical / ribosomal RNA / tRNA",
  (d) (per mod-type file) have no modifications of that mod type.
- No Markdown output, only TSV.
"""

import argparse
import pandas as pd
import re
import sys
import os


VALID_DTYPES = ("4mC", "5mC", "6mA")


def normalize_dtype(dtype: str) -> str:
    lookup = {
        "4mc": "4mC",
        "5mc": "5mC",
        "6ma": "6mA",
    }
    key = str(dtype).strip().lower()
    if key not in lookup:
        raise ValueError(
            f"Unsupported DATATYPE '{dtype}'. Valid values are: {', '.join(VALID_DTYPES)}"
        )
    return lookup[key]


def parse_args():
    p = argparse.ArgumentParser(
        description="Map mods to regions, annotate with product and stats, using +/– for diffs."
    )
    p.add_argument(
        "-i", dest="inputs", nargs=3, action="append", metavar=("FILE", "TREATMENT", "DATATYPE"),
        required=True,
        help="Mod file (BED 0-based), treatment name, and modification type (4mC, 5mC, or 6mA)"
    )
    p.add_argument(
        "-r", dest="regions", nargs=2, action="append", metavar=("FILE", "TYPE"),
        required=True,
        help="Region file and type: 'gene' (GFF3), 'upstream' (BED) or 'terminator' (BED)"
    )
    p.add_argument(
        "-d", dest="diffstats", nargs=2, action="append", default=[], metavar=("FILE", "GROUP"),
        help="CSV of differential stats (with locus_tag,log2FoldChange,…) and group name"
    )
    p.add_argument(
        "--gff", dest="gff3", required=True,
        help="Master GFF3 with product= annotations keyed by locus_tag"
    )
    p.add_argument(
        "-u", "--operon-map", dest="operon_map", default=None,
        help="TSV (no header): operon_id<TAB>locus_tag. Used to annotate upstream/terminator "
             "functions by listing all operon gene products and to compute per-group Log2FC ranges."
    )
    p.add_argument(
        "-o", dest="output", default=None,
        help=("Output TSV base name. Writes a combined raw table here, and also writes "
              "<base>_4mC_filtered.tsv, <base>_5mC_filtered.tsv, and <base>_6mA_filtered.tsv. "
              "If omitted, only the combined raw table is printed to stdout.")
    )
    return p.parse_args()


def load_operon_map(path: str) -> dict:
    """
    Read two-column TSV (no header): operon_id<TAB>locus_tag
    Returns: dict[operon_id] -> list[locus_tag] (order preserved).
    """
    df = pd.read_csv(path, sep="\t", header=None, names=["operon_id", "locus_tag"], dtype=str)
    df = df.dropna(subset=["operon_id", "locus_tag"])
    groups = df.groupby("operon_id")["locus_tag"].apply(list)
    return groups.to_dict()


def extract_locus_tag(attrs: str) -> str:
    m = re.search(r"locus_tag=([^;]+)", attrs)
    return m.group(1) if m else None


def extract_name(attrs: str) -> str:
    m = re.search(r"Name=([^;]+)", attrs)
    return m.group(1) if m else None


def parse_gff_products(path: str) -> dict:
    entries = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            lt = extract_locus_tag(parts[8])
            prod_m = re.search(r"product=([^;]+)", parts[8])
            prod = prod_m.group(1) if prod_m else None
            entries.append((lt, prod))
    prod_map = {}
    for i, (lt, prod) in enumerate(entries):
        if not lt or lt in prod_map:
            continue
        if prod:
            prod_map[lt] = prod
        else:
            if i + 1 < len(entries) and entries[i + 1][0] == lt and entries[i + 1][1]:
                prod_map[lt] = entries[i + 1][1]
    return prod_map


def load_regions(path: str, rtype: str) -> pd.DataFrame:
    rt = rtype.lower().rstrip("s")
    if rt == "gene":
        cols = ['contig', 'source', 'feat', 'start', 'end', 'score', 'strand', 'phase', 'attrs']
        df = pd.read_csv(path, sep="\t", comment="#", names=cols,
                         header=None, dtype=str)
        df['locus_tag'] = df['attrs'].apply(extract_locus_tag)
        df['name'] = df['attrs'].apply(extract_name)
        mask_gene = df['feat'].isin(["gene", "pseudogene"])
        mask_rna = df['source'] == 'RNA'
        df = df[mask_gene | mask_rna]
        df['locus_tag'] = df['locus_tag'].fillna(df['name'])
        df = df.dropna(subset=['locus_tag'])
        df['start'] = df['start'].astype(int) - 1
        df['end'] = df['end'].astype(int)
        df['region'] = df.apply(
            lambda x: f"gene:{x['locus_tag']}" if x['source'] == 'RNA'
            else f"{rt};{x['locus_tag']}",
            axis=1
        )
        df['region_type'] = rt
        return df[['contig', 'start', 'end', 'strand', 'region', 'locus_tag', 'region_type']]

    elif rt in ("upstream", "terminator"):
        cols = ['contig', 'start', 'end', 'locus_tag', 'score', 'strand']
        df = pd.read_csv(path, sep="\t", comment="#", names=cols,
                         header=None, dtype=str)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df['region'] = rt + ";" + df['locus_tag']
        df['region_type'] = rt
        return df[['contig', 'start', 'end', 'strand', 'region', 'locus_tag', 'region_type']]

    else:
        raise ValueError(f"Unknown region type: {rtype!r}")

def load_mods(path: str) -> pd.DataFrame:
    """
    Load a per-base DIFF table and return contig/start/end/strand/diff.
    Supports:
      (A) New/filtered files WITH header:
          columns include: chrom,start,end,strand,control_pct_modified,
          treatment_pct_modified,differential_pct_modification
      (B) Legacy 15-col files WITHOUT header
          (control at col index 12, treatment at col index 13, diff at col index 14)

    'diff' is returned as:
        differential_pct_modification[control_pct_modified,treatment_pct_modified]
    """
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)

    if "differential_pct_modification" in df.columns:
        rename_map = {
            "chrom": "contig",
            "differential_pct_modification": "diff",
        }
        if "strand" not in df.columns and df.shape[1] >= 6:
            df["strand"] = df.iloc[:, 5]

        df = df.rename(columns=rename_map)
        needed = [
            "contig",
            "start",
            "end",
            "strand",
            "diff",
            "control_pct_modified",
            "treatment_pct_modified",
        ]
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise ValueError(f"{path}: missing required columns {missing} in headered file")

    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, low_memory=False)
        if df.shape[1] < 15:
            raise ValueError(f"{path}: need ≥15 columns, got {df.shape[1]}")
        df = df.rename(columns={
            0: "contig",
            1: "start",
            2: "end",
            5: "strand",
            12: "control_pct_modified",
            13: "treatment_pct_modified",
            14: "diff",
        })

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["diff"] = pd.to_numeric(df["diff"], errors="coerce")

    df = df.dropna(subset=[
        "contig",
        "start",
        "end",
        "strand",
        "diff",
        "control_pct_modified",
        "treatment_pct_modified",
    ])

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["strand"] = df["strand"].astype(str)

    df["diff"] = df.apply(
        lambda x: (
            f"{x['diff']}["
            f"{str(x['control_pct_modified']).strip()},"
            f"{str(x['treatment_pct_modified']).strip()}]"
        ),
        axis=1,
    )

    return df[["contig", "start", "end", "strand", "diff"]]

def format_diff(d: float) -> str:
    s = str(d)
    return s if s.startswith('-') else f"+{s}"


def format_fc(d: float) -> str:
    return f"{d:+g}"


def is_bad_function_string(s: str) -> bool:
    """
    Return True if a function string should be considered 'bad' for filtering:
    - contains 'hypothetical'
    - or 'ribosomal rna'
    - or 'trna'
    Case-insensitive; robust against prefixes like '5S ribosomal RNA',
    but does NOT filter things like 'ribosomal protection-like ABC-F family protein'
    because that lacks 'RNA'.
    """
    if s is None:
        return False
    s_l = s.lower()
    if "hypothetical" in s_l:
        return True
    if "ribosomal rna" in s_l:
        return True
    if "trna" in s_l:
        return True
    return False


def filter_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply filtering rules:

    (a) Drop rows that don't have ALL Log2FC values (any empty "" or NaN).
    (b) If region starts with 'gene;' and function is hypothetical / ribosomal RNA / tRNA -> drop.
    (c) If region starts with 'upstream;' and function list *all* bad (as in (b)) -> drop.
        - If function is empty / whitespace / only empty entries -> keep.
        - If at least one entry is NOT bad -> keep.

    NOTE: Rule (d) (must have at least one mod of the relevant type) is applied
    outside this function, per mod-type subset.
    """
    if df is None or df.empty:
        return df

    log_cols = [c for c in df.columns if c.endswith("-Log2FC")]

    if log_cols:
        tmp = df[log_cols].fillna("")
        mask_has_all = (tmp != "").all(axis=1)
    else:
        mask_has_all = pd.Series(True, index=df.index)

    func = df.get("function", pd.Series("", index=df.index)).fillna("")

    is_gene_region = df["region"].str.startswith("gene;")
    bad_gene = is_gene_region & func.apply(is_bad_function_string)

    def upstream_is_bad(row) -> bool:
        region = row["region"]
        if not isinstance(region, str) or not region.startswith("upstream;"):
            return False

        f = row.get("function", "")
        if f is None:
            f = ""
        f = str(f)

        if f.strip() == "":
            return False

        entries = [x.strip() for x in f.split(";") if x.strip() != ""]
        if not entries:
            return False

        return all(is_bad_function_string(e) for e in entries)

    bad_upstream = df.apply(upstream_is_bad, axis=1)

    mask_keep = mask_has_all & (~bad_gene) & (~bad_upstream)

    return df[mask_keep].copy()

def print_region_overlap_summary(dfm: pd.DataFrame, mods: dict, labels_sorted):
    """
    Print a small terminal table showing, for each mod file (row),
    the % of input significant bases that overlap each region type,
    plus the % overlapping any provided region at all.

    A base is counted once per region type, even if it overlaps multiple
    regions of that type. In 'all_regions', a base is counted once total.
    """
    base_cols = ["contig", "start", "end", "strand"]

    total_by_label = {
        label: mods[label][base_cols].drop_duplicates().shape[0]
        for label in labels_sorted
    }

    hits = dfm[["treatment", "region_type"] + base_cols].drop_duplicates()

    rows = []
    for label in labels_sorted:
        total = total_by_label.get(label, 0)
        sub = hits[hits["treatment"] == label]

        row = {
            "mod": label,
            "n_bases": total,
        }

        for rt in ("gene", "upstream", "terminator"):
            n_rt = sub[sub["region_type"] == rt][base_cols].drop_duplicates().shape[0]
            row[rt] = 100.0 * n_rt / total if total else 0.0

        n_any = sub[base_cols].drop_duplicates().shape[0]
        row["all_regions"] = 100.0 * n_any / total if total else 0.0

        rows.append(row)

    summary = pd.DataFrame(rows)

    for col in ["gene", "upstream", "terminator", "all_regions"]:
        summary[col] = summary[col].map(lambda x: f"{x:.1f}%")

    print("\nRegion-overlap summary (% of input significant bases)", file=sys.stderr)
    print(summary.to_string(index=False), file=sys.stderr)


def main():
    args = parse_args()

    region_dfs = [load_regions(f, t) for f, t in args.regions]
    regions = pd.concat(region_dfs, ignore_index=True)

    product_map = parse_gff_products(args.gff3)

    operon_map = load_operon_map(args.operon_map) if args.operon_map else {}

    mods, labels = {}, []
    label_info = {}
    treatment_order = []

    for path, treat, dtype_raw in args.inputs:
        dtype = normalize_dtype(dtype_raw)
        label = f"{treat}_{dtype}"
        labels.append(label)
        mods[label] = load_mods(path)

        if treat not in treatment_order:
            treatment_order.append(treat)
        label_info[label] = {"treatment": treat, "dtype": dtype}

    recs = []
    for label, dfm in mods.items():
        for _, m in dfm.iterrows():
            same_region = (regions.contig == m.contig)
            overlap = (regions.start <= m.end) & (regions.end >= m.start)
            for _, r in regions[same_region & overlap].iterrows():
                recs.append({
                    'region': r.region,
                    'locus_tag': r.locus_tag,
                    'region_type': r.region_type,
                    'treatment': label,
                    'contig': m.contig,
                    'start': m.start,
                    'end': m.end,
                    'strand': m.strand,
                    'region_strand': r.strand,
                    'diff': m['diff']
                })

    if not recs:
        sys.stderr.write("No modifications fell within any region.\n")
        sys.exit(1)

    dfm = pd.DataFrame(recs)

    dfm['shorthand'] = dfm.apply(
        lambda x: f"{x['contig']}:{x['start']}-{x['end']}({x['strand']}):{format_diff(x['diff'])}",
        axis=1
    )

    def combine(group):
        same = group[group['strand'] == group['region_strand']]['shorthand']
        opp = group[group['strand'] != group['region_strand']]['shorthand']
        s_same = ";".join(same)
        s_opp = ";".join(opp)
        if s_same and s_opp:
            return f"{s_same}|||{s_opp}"
        elif s_same:
            return s_same
        elif s_opp:
            return f"|||{s_opp}"
        else:
            return ""

    pivot = (
        dfm
        .groupby(['region', 'locus_tag', 'treatment'])
        .apply(combine)
        .unstack('treatment')
        .fillna('')
        .reset_index()
    )

    pivot.insert(
        1,
        'function',
        pivot['locus_tag'].map(lambda lt: product_map.get(lt, ""))
    )
    mask_gene_rna = pivot['region'].str.startswith('gene:')
    pivot.loc[mask_gene_rna, 'function'] = pivot.loc[mask_gene_rna, 'locus_tag']

    if operon_map:
        mask_ut = pivot['region'].str.startswith(('upstream;', 'terminator;'))

        def operon_functions(row):
            opid = row['locus_tag']
            loci = operon_map.get(opid, [])
            if not loci:
                return row['function']
            prods = []
            for lt in loci:
                prod = product_map.get(lt, "")
                prods.append(prod if prod else lt)
            return ";".join(prods)

        pivot.loc[mask_ut, 'function'] = pivot.loc[mask_ut].apply(operon_functions, axis=1)

    coord_df = regions[['region', 'contig', 'start', 'end', 'strand']].drop_duplicates()
    pivot = pivot.merge(coord_df, on='region', how='left')
    pivot.insert(
        2,
        'coordinates',
        pivot.apply(lambda x: f"{x['contig']}:{x['start']}-{x['end']}({x['strand']})", axis=1)
    )
    pivot = pivot.drop(columns=['contig', 'start', 'end', 'strand'])

    dtype_priority = {"4mC": 0, "5mC": 1, "6mA": 2}

    def label_sort_key(lab: str):
        info = label_info.get(lab, {})
        treat = info.get("treatment", "")
        dtype = info.get("dtype", "")
        t_idx = treatment_order.index(treat) if treat in treatment_order else len(treatment_order)
        d_idx = dtype_priority.get(dtype, 99)
        return (t_idx, d_idx, lab)

    labels_sorted = sorted(labels, key=label_sort_key)

    out = pivot[['region', 'function', 'coordinates', 'locus_tag'] + labels_sorted].copy()

    for stats_file, group in args.diffstats:
        stats_df = pd.read_csv(stats_file)
        if 'locus_tag' not in stats_df.columns:
            stats_df = stats_df.rename(columns={stats_df.columns[0]: 'locus_tag'})
        if 'log2FoldChange' not in stats_df.columns:
            raise ValueError(f"{stats_file}: missing 'log2FoldChange'")

        series = pd.to_numeric(
            stats_df.set_index('locus_tag')['log2FoldChange'],
            errors='coerce'
        ).to_dict()

        col = f"{group}-Log2FC"

        def map_log2fc(row):
            reg = row['region']
            if operon_map and (reg.startswith('upstream;') or reg.startswith('terminator;')):
                opid = row['locus_tag']
                loci = operon_map.get(opid, [])
                vals = [series.get(lt) for lt in loci]
                vals = [v for v in vals if v is not None and pd.notnull(v)]
                if not vals:
                    return ""
                mn, mx = min(vals), max(vals)
                return f"{format_fc(mn)} to {format_fc(mx)}"
            v = series.get(row['locus_tag'])
            return format_fc(v) if v is not None and pd.notnull(v) else ""

        out[col] = out.apply(map_log2fc, axis=1)

    out = out.drop(columns=['locus_tag'])

    base_cols = ['region', 'function', 'coordinates']
    log_cols = [c for c in out.columns if c.endswith("-Log2FC")]

    def log_sort_key(col: str):
        group = col[:-7]
        if group in treatment_order:
            return (0, treatment_order.index(group))
        else:
            return (1, group)

    log_cols_sorted = sorted(log_cols, key=log_sort_key)

    all_cols = base_cols + labels_sorted + log_cols_sorted
    out_all = out[all_cols].copy()
    print_region_overlap_summary(dfm, mods, labels_sorted)
    
    if not args.output:
        out_all.to_csv(sys.stdout, sep="\t", index=False)
        return

    base, _ = os.path.splitext(args.output)

    out_all.to_csv(args.output, sep="\t", index=False)

    mod_cols_4mc = []
    mod_cols_5mc = []
    mod_cols_6ma = []
    for lab in labels_sorted:
        info = label_info.get(lab, {})
        dtype = info.get("dtype")
        if dtype == "4mC" and lab in out_all.columns:
            mod_cols_4mc.append(lab)
        elif dtype == "5mC" and lab in out_all.columns:
            mod_cols_5mc.append(lab)
        elif dtype == "6mA" and lab in out_all.columns:
            mod_cols_6ma.append(lab)

    def build_filtered(out_all: pd.DataFrame, mod_cols, suffix: str):
        if not mod_cols:
            return

        df = out_all[base_cols + mod_cols + log_cols_sorted].copy()

        mods_only = df[mod_cols].fillna("")
        has_any_mod = (mods_only != "").any(axis=1)
        df = df[has_any_mod].copy()
        if df.empty:
            return

        df_filt = filter_table(df)
        if df_filt.empty:
            return

        df_filt.to_csv(f"{base}{suffix}", sep="\t", index=False)

    build_filtered(out_all, mod_cols_4mc, "_4mC_filtered.tsv")
    build_filtered(out_all, mod_cols_5mc, "_5mC_filtered.tsv")
    build_filtered(out_all, mod_cols_6ma, "_6mA_filtered.tsv")


if __name__ == "__main__":
    main()