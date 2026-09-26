#!/usr/bin/env python3
"""
3D Bivariate Histogram Plot Script

This script generates a 3D bivariate histogram for each gene (locus tag)
based on BLAST output, metadata, and an index of locus tags provided
as a FASTA file (one header per locus tag).

Changes from the original:
  1) --index now points to a FASTA file; each header is a locus tag to plot.
  2) -o/--output is an output FOLDER; the script writes one PNG per locus,
     named <locus_tag>.png, and nothing else.
  3) Per-plot titles are just the locus tag (no unique colors/descriptions).

Everything else remains as before.
"""

import argparse
import csv
from collections import defaultdict
import math
import os
import re
import matplotlib
matplotlib.use('Agg')  # For headless systems.
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # Enable 3D plotting.
import numpy as np

#############################
# Configuration Parameters  #
#############################

# General Data Filter Settings:
GENERAL_SETTINGS = {
    'start_year': 1950  # All data with a year before this value will be ignored.
}

# Wireframe (surface) settings:
WIREFRAME_SETTINGS = {
    'cmap': 'viridis',     # Colormap for the surface.
    'linewidth': 0.1,      # Thickness of the wireframe lines.
    'edgecolor': 'k',      # Color of the wireframe edges.
    'antialiased': True,   # Smoothing flag.
    'alpha': 0.88          # Transparency.
}

# Scatter marker (dot) settings:
SCATTER_SETTINGS = {
    's': 0.00005,          # Marker area.
    'c': 'k',              # Marker color (fill).
    'alpha': 0.8,          # Transparency.
    'marker': 'o',         # Marker shape (circle).
    'edgecolors': 'k'      # Edge color; same as fill for a fully filled circle.
}

# Trend line settings:
TREND_LINE_SETTINGS = {
    'color': 'r',          # Trend line color.
    'linewidth': 1.5       # Trend line thickness.
}
TREND_LINE_EPSILON = 0.5

# Axis settings (for X and Z axes):
AXIS_SETTINGS = {
    'x_label': 'Year',
    'z_label': 'Prevalence (%)'
}

# Identity (Y) axis settings:
IDENTITY_SETTINGS = {
    'range': (90, 100),              # Identity axis range (min, max) in percent.
    'ticks': np.arange(90, 101, 5),  # Tick marks for identity axis.
    'bins': 50                       # Number of bins (boxes) across the identity range.
}

# Subplot layout (kept for consistency; only single-plot saving is used now)
SUBPLOT_ADJUSTMENTS = {
    'left': 0.1,
    'right': 0.9,
    'top': 0.9,
    'bottom': 0.1,
    'wspace': 0.25,
    'hspace': 0.1
}

# ----------------------------------------------------------------

#############################
# Data Loading Functions    #
#############################

def load_metadata(meta_file):
    """
    Loads metadata from a TSV file with columns including at least:
      id, year

    Returns:
      meta: dict mapping accession (int) -> float(year)
      year_assemblies: dict mapping year (int) -> set(accession)
    """
    meta = {}
    year_assemblies = defaultdict(set)
    try:
        with open(meta_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    year_str = row.get("year", "").strip()
                    if not year_str or year_str.lower() == "unknown":
                        continue
                    year_float = float(year_str)
                    if year_float < GENERAL_SETTINGS['start_year']:
                        continue
                    id_str = row.get("id", "").strip()
                    if not id_str:
                        continue
                    accession = int(id_str)
                    if accession in meta:
                        if year_float < meta[accession]:
                            old_year = int(meta[accession])
                            year_assemblies[old_year].discard(accession)
                            meta[accession] = year_float
                            year_assemblies[int(year_float)].add(accession)
                    else:
                        meta[accession] = year_float
                        year_assemblies[int(year_float)].add(accession)
                except Exception:
                    continue
    except Exception as e:
        print(f"Error reading metadata file '{meta_file}': {e}")
    return meta, year_assemblies

def load_index_fasta(index_fa_path):
    """
    Loads locus tags from a FASTA file: each header line (>locus_tag ...)
    contributes one locus tag. Order is preserved, duplicates removed.
    """
    loci_ordered = []
    seen = set()
    try:
        with open(index_fa_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    tag = line[1:].strip().split()[0]
                    if tag and tag not in seen:
                        seen.add(tag)
                        loci_ordered.append(tag)
    except Exception as e:
        print(f"Error reading index FASTA '{index_fa_path}': {e}")
    return loci_ordered, set(loci_ordered)

def parse_blast_file(file_path, meta, locus_set):
    """
    Reads a BLAST output file with at least three tab-separated columns:
      - Query gene/locus in col 1 (e.g. "pgaptmp_001076::chromosome:..." or "pgaptmp_001076")
      - Subject (col 2): sseqid (e.g. "1234;contig_7" or "1234_7" or "lcl|1234_7" or "1234-7" or "1234")
      - % Identity (col 3)
    Filters to rows where the query locus is in locus_set and the subject accession
    can be mapped to an integer present in metadata.

    Returns list of tuples: (locus_tag, acc, year_float, pid)
    """
    hits = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 3:
                    continue

                gene_field = fields[0].strip()
                locus_tag = gene_field.split("::")[0] if "::" in gene_field else gene_field
                if locus_tag not in locus_set:
                    continue

                acc_field = fields[1].strip()
                try:
                    pid = float(fields[2])
                except ValueError:
                    continue

                # Minimal, format-tolerant extraction of integer accession from sseqid:
                # try common delimiters first; otherwise fall back to first integer in the string.
                acc_token = acc_field
                if ';' in acc_field:
                    acc_token = acc_field.split(';')[0]
                elif '_' in acc_field:
                    acc_token = acc_field.split('_')[0]
                elif '-' in acc_field:
                    acc_token = acc_field.split('-')[0]
                elif '|' in acc_field:
                    # e.g., lcl|1234_7 -> take the last piece
                    acc_token = acc_field.split('|')[-1]

                m = re.match(r'^(\d+)$', acc_token)
                if not m:
                    m = re.search(r'(\d+)', acc_field)
                if not m:
                    continue
                try:
                    acc = int(m.group(1))
                except Exception:
                    continue

                if acc in meta:
                    year_float = meta[acc]
                    hits.append((locus_tag, acc, year_float, pid))
    except Exception as e:
        print(f"Error reading BLAST file '{file_path}': {e}")
    return hits

def compute_xticks(start_year, end_year):
    """
    Computes x-axis tick positions starting at the nearest multiple of 10,
    then every 20 years until less than end_year.
    """
    lower_tick = math.ceil(start_year / 10) * 10
    ticks = []
    t = lower_tick
    while t < end_year:
        ticks.append(t)
        t += 20
    return ticks

#############################
# Single Gene Plot Function #
#############################

def plot_single_gene(gene, years_list, xticks, start_year, end_year,
                     bin_edges, bin_centers, meta, year_assemblies,
                     gene_year_assemblies_final, gene_identity_unique,
                     identity_max, identity_settings, ax=None):
    """
    Plots the 3D bivariate histogram for a single gene.
    If ax is None, a new figure and axis are created.
    Returns the axis containing the plot.
    """
    if ax is None:
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')
    ax.grid(False)
    
    Z = np.zeros((len(years_list), len(bin_edges)-1), dtype=float)
    data_found = False
    for i, yr in enumerate(years_list):
        total_in_year = len(year_assemblies.get(yr, set()))
        if total_in_year < 1:
            continue
        if yr not in gene_identity_unique[gene]:
            continue
        gene_count = len(gene_year_assemblies_final[gene][yr])
        prevalence = (gene_count / total_in_year) * 100
        identities = list(gene_identity_unique[gene][yr].values())
        if not identities:
            continue
        data_found = True
        if len(identities) == 1:
            val = identities[0]
            if val < bin_edges[0]:
                val = bin_edges[0]
            bin_idx = np.searchsorted(bin_edges, val) - 1
            bin_idx = max(bin_idx, 0)
            bin_idx = min(bin_idx, len(bin_edges)-2)
            Z[i, bin_idx] = prevalence
        else:
            clamped = [max(x, bin_edges[0]) for x in identities]
            hist, _ = np.histogram(clamped, bins=bin_edges)
            fractions = hist / len(clamped)
            scaled = fractions * prevalence
            Z[i, :] = scaled

    years_grid, id_grid = np.meshgrid(years_list, bin_centers, indexing='ij')
    if data_found and np.sum(Z) > 0:
        ax.plot_surface(
            years_grid,
            id_grid,
            Z,
            cmap=WIREFRAME_SETTINGS['cmap'],
            linewidth=WIREFRAME_SETTINGS['linewidth'],
            edgecolor=WIREFRAME_SETTINGS['edgecolor'],
            antialiased=WIREFRAME_SETTINGS['antialiased'],
            alpha=WIREFRAME_SETTINGS['alpha']
        )
        ax.scatter(
            years_grid.ravel(),
            id_grid.ravel(),
            Z.ravel(),
            s=SCATTER_SETTINGS['s'],
            c=SCATTER_SETTINGS['c'],
            alpha=SCATTER_SETTINGS['alpha'],
            marker=SCATTER_SETTINGS['marker'],
            edgecolors=SCATTER_SETTINGS['edgecolors']
        )
        # Trend line over 3-year intervals.
        trend_data = {}
        for yr in range(start_year, end_year + 1, 3):
            assemblies = year_assemblies.get(yr, set())
            if len(assemblies) < 1:
                continue
            if yr in gene_year_assemblies_final[gene]:
                gene_count = len(gene_year_assemblies_final[gene][yr])
                prevalence_yr = (gene_count / len(assemblies)) * 100
                trend_data[yr] = prevalence_yr
        if trend_data:
            first_valid = min(trend_data.keys())
            trend_x = []
            trend_z = []
            last_valid_value = None
            for yr in range(first_valid, end_year + 1, 3):
                if yr in trend_data:
                    last_valid_value = trend_data[yr]
                if last_valid_value is not None:
                    trend_x.append(yr)
                    trend_z.append(last_valid_value)
            if len(trend_x) >= 2:
                trend_y = [identity_max - 0.5] * len(trend_x)
                ax.plot(trend_x, trend_y, trend_z,
                        color=TREND_LINE_SETTINGS['color'],
                        linewidth=TREND_LINE_SETTINGS['linewidth'])
    else:
        ax.text2D(0.5, 0.5, "No data", transform=ax.transAxes,
                  horizontalalignment='center', verticalalignment='center')
    
    # Axis control annotations
    label_fs = 16
    tick_fs = 14

    ax.set_xlabel(AXIS_SETTINGS['x_label'], fontsize=label_fs, labelpad=10)
    ax.set_ylabel('Identity (%)', fontsize=label_fs, labelpad=10)
    ax.set_zlabel(AXIS_SETTINGS['z_label'], fontsize=label_fs, labelpad=12)
    
    ax.tick_params(axis='x', labelsize=tick_fs, rotation=0)
    ax.tick_params(axis='y', labelsize=tick_fs)
    ax.tick_params(axis='z', labelsize=tick_fs)

    ax.set_xticks(xticks)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(IDENTITY_SETTINGS['range'][0], IDENTITY_SETTINGS['range'][1])
    ax.set_yticks(IDENTITY_SETTINGS['ticks'])

    ax.set_zticks([0, 50, 100])
    
    return ax

#############################
# Main Plotting Function    #
#############################

def main():
    parser = argparse.ArgumentParser(
        description="3D bivariate histogram plots per locus tag."
    )
    parser.add_argument('-meta', '--metadata', required=True,
                        help='Path to metadata TSV file')
    parser.add_argument('-in', '--input', required=True,
                        help='Path to BLAST output file')
    parser.add_argument('--index', required=True,
                        help='Path to FASTA listing locus tags (one header per locus)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads (unused, for compatibility)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output folder in which one <locus_tag>.png is written per locus')
    args = parser.parse_args()
    
    # Prepare output folder
    out_dir = os.path.abspath(args.output)
    os.makedirs(out_dir, exist_ok=True)

    meta, year_assemblies = load_metadata(args.metadata)
    if not meta:
        print("Warning: No valid metadata loaded. Exiting.")
        return

    locus_list, locus_set = load_index_fasta(args.index)
    if not locus_list:
        print("Warning: No locus tags found in --index FASTA. Exiting.")
        return
    
    all_hits = parse_blast_file(args.input, meta, locus_set)
    if not all_hits:
        print("Warning: No valid BLAST hits found.")

    # Aggregate by locus/year and keep best identity per accession per year
    gene_year_assemblies_final = defaultdict(lambda: defaultdict(set))
    gene_identity_unique = defaultdict(lambda: defaultdict(dict))
    for locus_tag, acc, year_float, pid in all_hits:
        yr = int(year_float)
        gene_year_assemblies_final[locus_tag][yr].add(acc)
        if acc not in gene_identity_unique[locus_tag][yr]:
            gene_identity_unique[locus_tag][yr][acc] = pid
        else:
            gene_identity_unique[locus_tag][yr][acc] = max(gene_identity_unique[locus_tag][yr][acc], pid)
    
    try:
        computed_start_year = min(int(y) for y in meta.values())
        end_year = max(int(y) for y in meta.values())
    except ValueError:
        computed_start_year, end_year = 2000, 2025
    start_year = max(GENERAL_SETTINGS['start_year'], computed_start_year)
    
    years_list = list(range(start_year, end_year + 1))
    xticks = compute_xticks(start_year, end_year)
    
    forced_min_identity = IDENTITY_SETTINGS['range'][0]
    identity_max = IDENTITY_SETTINGS['range'][1]
    num_bins = IDENTITY_SETTINGS['bins']
    bin_edges = np.linspace(forced_min_identity, identity_max, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # --- Per-locus plots directly into output folder ---
    for gene in locus_list:
        fig_single = plt.figure(figsize=(7,7))
        ax_single = fig_single.add_subplot(111, projection='3d')
        ax_single = plot_single_gene(
            gene, years_list, xticks, start_year, end_year,
            bin_edges, bin_centers, meta, year_assemblies,
            gene_year_assemblies_final, gene_identity_unique,
            identity_max, IDENTITY_SETTINGS, ax=ax_single
        )
        # Title: locus tag only
        ax_single.set_title(gene, fontsize=10)
        gene_filename = os.path.join(out_dir, f"{gene}.png")
        try:
            fig_single.savefig(gene_filename, dpi=600)
        except Exception as e:
            print(f"Error saving plot for gene '{gene}': {e}")
        finally:
            plt.close(fig_single)

if __name__ == '__main__':
    main()
