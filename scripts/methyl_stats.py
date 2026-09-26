#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarise DIFF methylation tables across conditions.

Usage example:
python methyl_stats.py \
  -i path/4mC_5mC_VANC.filtered.bed path/4mC_5mC_VANC.DIFF.bed VANC_4mC_5mC \
  -i path/6mA_VANC.filtered.bed     path/6mA_VANC.DIFF.bed      VANC_6mA \
  -o methyl_report.tsv
"""

import argparse
import csv
import math
from statistics import mean, median

def safe_int(x, default=None):
    try:
        return int(x)
    except Exception:
        return default


def load_bed_regions(path):
    """
    Load a BED file with at least 3 columns: chrom, start, end.
    Strand is ignored. Returns: dict {chrom: [(start, end), ...]} with
    intervals sorted by start.
    """
    regions = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = safe_int(parts[1])
            end = safe_int(parts[2])
            if start is None or end is None:
                continue
            regions.setdefault(chrom, []).append((start, end))
    for chrom in regions:
        regions[chrom].sort(key=lambda x: x[0])
    return regions


def load_bed_regions_stranded(path):
    """
    Load a BED file with at least 3 columns: chrom, start, end, plus
    optional strand in column 6 (index 5). Returns:
        {chrom: {strand: [(start, end), ...]}}
    with intervals sorted by start for each strand.
    """
    regions = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = safe_int(parts[1])
            end = safe_int(parts[2])
            if start is None or end is None:
                continue
            strand = parts[5] if len(parts) > 5 else "."
            chrom_map = regions.setdefault(chrom, {})
            chrom_map.setdefault(strand, []).append((start, end))
    for chrom in regions:
        for strand in regions[chrom]:
            regions[chrom][strand].sort(key=lambda x: x[0])
    return regions


def pos_in_intervals(regions, chrom, pos):
    """
    Return True if position 'pos' on 'chrom' lies within any interval
    in 'regions' (BED-like, 0-based half-open: start <= pos < end).
    Uses binary search per chromosome. Strand-agnostic.
    """
    intervals = regions.get(chrom)
    if not intervals:
        return False
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        start, end = intervals[mid]
        if pos < start:
            hi = mid
        elif pos >= end:
            lo = mid + 1
        else:
            return True
    return False


def pos_in_intervals_stranded(regions, chrom, pos, strand):
    """
    Return True if position 'pos' on 'chrom' lies within any interval
    in 'regions' for the given strand. 'regions' is of the form:
        {chrom: {strand: [(start, end), ...]}}
    Uses binary search per chromosome+strand.
    """
    chrom_map = regions.get(chrom)
    if not chrom_map:
        return False
    intervals = chrom_map.get(strand)
    if not intervals:
        return False
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        start, end = intervals[mid]
        if pos < start:
            hi = mid
        elif pos >= end:
            lo = mid + 1
        else:
            return True
    return False



def read_diff_file(path):
    """
    Load a DIFF file (tab-delimited with header).
    Returns (rows, header). Each row is a list of strings (no casting here).
    Assumptions about columns (0-based indices):
      7  -> control_total
      9  -> treatment_total
      12 -> control_pct_modified  (fraction, e.g., 0.30)
      13 -> treatment_pct_modified (fraction)
      -1 -> diff (treatment - control), fraction
    """
    rows = []
    with open(path, newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd, None)
        for r in rd:
            if r:
                rows.append(r)
    return rows, header

def safe_float(x, default=None):
    try:
        v = float(x)
        return v if math.isfinite(v) else default
    except Exception:
        return default

def basic_stats(values):
    """
    Return (count, mean, median, min, max) for a numeric list.
    Empty-safe: returns (0, None, None, None, None) for empty lists.
    """
    if not values:
        return 0, None, None, None, None
    return len(values), mean(values), median(values), min(values), max(values)

def compute_metrics(
    sig_rows,
    all_rows,
    gene_regions=None,
    upstream_regions=None,
    gene_regions_stranded=None,
    upstream_regions_stranded=None,
):
    """
    Compute a compact, high-signal set of metrics.

    Inputs:
      sig_rows: rows from the 'significant' file (already filtered upstream)
      all_rows: rows from the 'all loci' DIFF file
      gene_regions: dict {chrom: [(start, end), ...]} for gene bodies (strand-agnostic)
      upstream_regions: same for upstream regions (strand-agnostic)
      gene_regions_stranded: {chrom: {strand: [(start, end), ...]}} for gene bodies
      upstream_regions_stranded: same for upstream regions

    Column indices (0-based):
      control_total           = 7
      treatment_total         = 9
      control_pct_modified    = 12
      treatment_pct_modified  = 13
      diff (treatment-control)= -1  (last column), fraction units

    Positional indices (assumed):
      chrom = 0
      start = 1   (single-base loci: start is the position)
      strand = 5  (BED-like, '+' or '-'; ignored if absent)
    """
    idx_ctot, idx_ttot = 7, 9
    idx_cpct, idx_tpct = 12, 13
    idx_diff = -1

    # ---------- SIGNIFICANT set ----------
    sig_diffs = []
    for r in sig_rows:
        d = safe_float(r[idx_diff])
        if d is not None:
            sig_diffs.append(d)

    abs_sig_diffs = [abs(x) for x in sig_diffs]
    _, mean_abs_sig, med_abs_sig, _, max_abs_sig = basic_stats(abs_sig_diffs)

    # ---------- ALL set ----------
    all_diffs = []
    all_ctot = []
    all_ttot = []
    zeros_in_control = 0

    thr20 = 0  # count(|diff| > 0.20)

    for r in all_rows:
        d = safe_float(r[idx_diff])
        if d is not None:
            all_diffs.append(d)
            if abs(d) > 0.20:
                thr20 += 1

        ctot = safe_float(r[idx_ctot])
        ttot = safe_float(r[idx_ttot])
        if ctot is not None:
            all_ctot.append(ctot)
        if ttot is not None:
            all_ttot.append(ttot)

        cpm = safe_float(r[idx_cpct])
        # "0% methylation in control" means exactly zero fraction
        if cpm is not None and cpm == 0.0:
            zeros_in_control += 1

    abs_all_diffs = [abs(x) for x in all_diffs]
    _, mean_abs_all, med_abs_all, _, _ = basic_stats(abs_all_diffs)

    # Counts
    n_sig = len(sig_rows)
    n_all = len(all_rows)

    # Percent helpers
    def pct(v):
        return None if v is None else f"{v:.2f}"

    def pct_from_fraction(v):
        return None if v is None else f"{(v*100.0):.2f}"

    # Percent with 0% in control among ALL loci
    pct_zero_control = (zeros_in_control / n_all * 100.0) if n_all else None

    # Percent of ALL loci with |diff|>20%
    pct_thr20 = (thr20 / n_all * 100.0) if n_all else None

    # Coverage proxies (means across ALL loci)
    mean_ctot = mean(all_ctot) if all_ctot else None
    mean_ttot = mean(all_ttot) if all_ttot else None

    # ---------- SIGNIFICANT set: localisation (strand-agnostic and strand-aware) ----------
    in_gene_any = 0
    in_upstream_any = 0
    in_gene_or_upstream_any = 0

    in_gene_any_strand = 0
    in_upstream_any_strand = 0
    in_gene_or_upstream_any_strand = 0

    if n_sig and (gene_regions is not None or upstream_regions is not None):
        for r in sig_rows:
            if not r:
                continue
            chrom = r[0]
            pos = safe_int(r[1])
            if pos is None:
                continue

            strand = r[5] if len(r) > 5 else None

            # Strand-agnostic overlaps
            in_gene = pos_in_intervals(gene_regions, chrom, pos) if gene_regions is not None else False
            in_upstream = pos_in_intervals(upstream_regions, chrom, pos) if upstream_regions is not None else False

            if in_gene:
                in_gene_any += 1
            if in_upstream:
                in_upstream_any += 1
            if in_gene or in_upstream:
                in_gene_or_upstream_any += 1

            # Strand-aware overlaps (only if strand and strand-specific regions provided)
            in_gene_s = False
            in_upstream_s = False
            if strand in ("+", "-"):
                if gene_regions_stranded is not None:
                    in_gene_s = pos_in_intervals_stranded(gene_regions_stranded, chrom, pos, strand)
                if upstream_regions_stranded is not None:
                    in_upstream_s = pos_in_intervals_stranded(upstream_regions_stranded, chrom, pos, strand)

            if in_gene_s:
                in_gene_any_strand += 1
            if in_upstream_s:
                in_upstream_any_strand += 1
            if in_gene_s or in_upstream_s:
                in_gene_or_upstream_any_strand += 1

        pct_gene_any = (in_gene_any / n_sig * 100.0) if n_sig else None
        pct_upstream_any = (in_upstream_any / n_sig * 100.0) if n_sig else None
        pct_gene_or_upstream_any = (in_gene_or_upstream_any / n_sig * 100.0) if n_sig else None

        pct_gene_any_strand = (in_gene_any_strand / n_sig * 100.0) if n_sig else None
        pct_upstream_any_strand = (in_upstream_any_strand / n_sig * 100.0) if n_sig else None
        pct_gene_or_upstream_any_strand = (in_gene_or_upstream_any_strand / n_sig * 100.0) if n_sig else None
    else:
        pct_gene_any = None
        pct_upstream_any = None
        pct_gene_or_upstream_any = None
        pct_gene_any_strand = None
        pct_upstream_any_strand = None
        pct_gene_or_upstream_any_strand = None

    metrics = {
        "NUM_significant_bases": n_sig,
        "NUM_all_bases": n_all,
        "Mean_abs_diff_%_(sig)": pct_from_fraction(mean_abs_sig),
        "Median_abs_diff_%_(sig)": pct_from_fraction(med_abs_sig),
        "Max_abs_diff_%_(sig)": pct_from_fraction(max_abs_sig),
        "Mean_abs_diff_%_(all)": pct_from_fraction(mean_abs_all),
        "Median_abs_diff_%_(all)": pct_from_fraction(med_abs_all),
        "Pct_bases_control_0%_(all)": pct(pct_zero_control),
        "Pct_|diff|>20%_(all)": pct(pct_thr20),
        "Mean_control_total_reads": f"{mean_ctot:.2f}" if mean_ctot is not None else None,
        "Mean_treatment_total_reads": f"{mean_ttot:.2f}" if mean_ttot is not None else None,
        # Strand-agnostic localisation
        "Pct_sig_bases_in_genes_any": pct(pct_gene_any),
        "Pct_sig_bases_in_upstream_any": pct(pct_upstream_any),
        "Pct_sig_bases_in_gene_or_upstream_any": pct(pct_gene_or_upstream_any),
        # Strand-aware localisation
        "Pct_sig_bases_in_genes_any_strand_match": pct(pct_gene_any_strand),
        "Pct_sig_bases_in_upstream_any_strand_match": pct(pct_upstream_any_strand),
        "Pct_sig_bases_in_gene_or_upstream_any_strand_match": pct(pct_gene_or_upstream_any_strand),
    }
    return metrics

def main():
    ap = argparse.ArgumentParser(
        description="Summarise DIFF methylation tables across conditions."
    )
    ap.add_argument(
        "-i", "--input", metavar=("SIGNIFICANT", "ALL", "LABEL"),
        nargs=3, action="append", required=True,
        help="Triplet: significant.tsv all.tsv label  (repeatable)"
    )
    ap.add_argument(
        "-o", "--out", required=True, help="Output TSV path."
    )
    ap.add_argument(
        "--gene", required=True, help="BED file with gene regions."
    )
    ap.add_argument(
        "--upstream", required=True, help="BED file with upstream regions."
    )
    args = ap.parse_args()

    # Load gene and upstream regions (strand-agnostic and strand-aware)
    gene_regions = load_bed_regions(args.gene)
    upstream_regions = load_bed_regions(args.upstream)
    gene_regions_stranded = load_bed_regions_stranded(args.gene)
    upstream_regions_stranded = load_bed_regions_stranded(args.upstream)

    labels = []
    per_label_metrics = {}

    # Compute metrics per label
    for sig_path, all_path, label in args.input:
        sig_rows, _ = read_diff_file(sig_path)
        all_rows, _ = read_diff_file(all_path)
        metrics = compute_metrics(
            sig_rows,
            all_rows,
            gene_regions=gene_regions,
            upstream_regions=upstream_regions,
            gene_regions_stranded=gene_regions_stranded,
            upstream_regions_stranded=upstream_regions_stranded,
        )
        labels.append(label)
        per_label_metrics[label] = metrics

    # Define output row order (clean & compact)
    metric_names = [
        "NUM_significant_bases",
        "NUM_all_bases",
        "Mean_abs_diff_%_(sig)",
        "Median_abs_diff_%_(sig)",
        "Max_abs_diff_%_(sig)",
        "Mean_abs_diff_%_(all)",
        "Median_abs_diff_%_(all)",
        "Pct_bases_control_0%_(all)",
        "Pct_|diff|>20%_(all)",
        "Mean_control_total_reads",
        "Mean_treatment_total_reads",
        "Pct_sig_bases_in_genes_any",
        "Pct_sig_bases_in_upstream_any",
        "Pct_sig_bases_in_gene_or_upstream_any",
        "Pct_sig_bases_in_genes_any_strand_match",
        "Pct_sig_bases_in_upstream_any_strand_match",
        "Pct_sig_bases_in_gene_or_upstream_any_strand_match",
    ]

    with open(args.out, "w", newline='') as fh:
        wr = csv.writer(fh, delimiter='\t')

        # Header
        wr.writerow(["Metric"] + labels)

        # Metrics table
        for metric in metric_names:
            row = [metric]
            for lab in labels:
                val = per_label_metrics[lab].get(metric)
                row.append(val if val is not None else "")
            wr.writerow(row)

        # Definitions (≤30 words each)
        wr.writerow([])
        wr.writerow(["# DEFINITIONS (≤30 words each)"])
        defs = {
            "NUM_significant_bases": "Count of loci in the pre-filtered ‘significant’ file (header removed).",
            "NUM_all_bases": "Total loci in the unfiltered genome-wide DIFF file (header removed).",
            "Mean_abs_diff_%_(sig)": "Average magnitude of methylation difference at significant loci (|treatment−control|), as percentage.",
            "Median_abs_diff_%_(sig)": "Median magnitude of methylation difference at significant loci, as percentage.",
            "Max_abs_diff_%_(sig)": "Largest magnitude methylation difference among significant loci, as percentage.",
            "Mean_abs_diff_%_(all)": "Average magnitude of methylation difference across all loci, as percentage.",
            "Median_abs_diff_%_(all)": "Median magnitude of methylation difference across all loci, as percentage.",
            "Pct_bases_control_0%_(all)": "Percentage of all loci with zero methylation in control.",
            "Pct_|diff|>20%_(all)": "Percentage of all loci with absolute methylation difference exceeding 20 percentage points.",
            "Mean_control_total_reads": "Mean control read depth (per-locus total reads) across all loci.",
            "Mean_treatment_total_reads": "Mean treatment read depth (per-locus total reads) across all loci.",
            "Pct_sig_bases_in_genes_any": "Percentage of significant loci that fall within any annotated gene body (strand-agnostic).",
            "Pct_sig_bases_in_upstream_any": "Percentage of significant loci that fall within any defined upstream region (strand-agnostic).",
            "Pct_sig_bases_in_gene_or_upstream_any": "Percentage of significant loci in either a gene body or an upstream region (strand-agnostic).",
            "Pct_sig_bases_in_genes_any_strand_match": "Percentage of significant loci in gene bodies where locus strand matches gene strand.",
            "Pct_sig_bases_in_upstream_any_strand_match": "Percentage of significant loci in upstream regions where locus strand matches region strand.",
            "Pct_sig_bases_in_gene_or_upstream_any_strand_match": "Percentage of significant loci in gene or upstream with matching strand.",
        }
        for k in metric_names:
            if k in defs:
                wr.writerow([k, defs[k]])


if __name__ == "__main__":
    main()
