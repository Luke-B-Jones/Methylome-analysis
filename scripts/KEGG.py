#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import math
import itertools
import matplotlib.ticker as ticker

from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerPatch

# =====================
#  Colors for each treatment: (up_color, down_color)
#  The legend will show a single "split box" for each treatment.
# =====================
TREATMENT_COLORS = {
    'H2O2': ('#1b9e77', '#8dcfbb'),
    'VANC': ('#7570b3', '#bab8d9'),
    'AMPI': ('#d95f02', '#ec8780'),
}

# Slightly increase the base font size
plt.rcParams.update({'font.size': 11})

class SplitColorPatch(Patch):
    """
    A dummy patch that holds two colors (up/down).
    We'll use a custom legend handler to draw it as a split box.
    """
    def __init__(self, up_color, down_color, label=None):
        super().__init__()
        self.up_color = up_color
        self.down_color = down_color
        self._label = label

class HandlerSplitColorPatch(HandlerPatch):
    """
    Custom legend handler that draws a single box split vertically
    into up_color (left) and down_color (right).
    """
    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        left_rect = Rectangle([xdescent, ydescent], width/2, height,
                              facecolor=orig_handle.up_color, transform=trans)
        right_rect = Rectangle([xdescent + width/2, ydescent], width/2, height,
                               facecolor=orig_handle.down_color, transform=trans)
        return [left_rect, right_rect]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Produce (1) a single diverging barplot (log scale) "
                    "and (2) a subplot bubble chart for KEGG pathways."
    )
    parser.add_argument('--data', nargs=2, action='append', required=True,
                        help='TSV file and treatment name (e.g. --data file.tsv H2O2). '
                             'File must have columns: treatment, KEGG_ID, SimpleDesc, direction, count.')
    parser.add_argument('-o', required=True,
                        help='Output base filename (e.g. kegg_plot.png). '
                             'A second figure is saved as *_subplots.png.')
    return parser.parse_args()

def read_and_aggregate(args):
    """ Read & aggregate data from all input files. """
    all_dfs = []
    for file, group in args.data:
        if not os.path.exists(file):
            sys.exit(f"Error: Data file '{file}' does not exist.")
        df = pd.read_csv(file, sep="\t")

        # Accept GSEA-style inputs:
        #   ID, Description and/or SimpleDesc, NES, core_enrichment (and/or setSize)
        # Normalize them to: treatment, KEGG_ID, SimpleDesc, direction, count
        if 'ID' not in df.columns:
            sys.exit("Error: Data file must include column 'ID' (KEGG pathway ID).")
        if 'NES' not in df.columns:
            sys.exit("Error: Data file must include column 'NES' to derive direction.")
        if ('core_enrichment' not in df.columns) and ('setSize' not in df.columns):
            sys.exit("Error: Need either 'core_enrichment' or 'setSize' to derive counts.")

        # KEGG_ID
        df = df.rename(columns={'ID': 'KEGG_ID'})
        df['KEGG_ID'] = df['KEGG_ID'].astype(str)

        # SimpleDesc: prefer provided SimpleDesc, else fall back to Description, else empty
        if 'SimpleDesc' in df.columns:
            simple_desc = df['SimpleDesc'].fillna('')
        elif 'Description' in df.columns:
            simple_desc = df['Description'].fillna('')
        else:
            simple_desc = ''

        # direction from NES (coerce to numeric first)
        nes = pd.to_numeric(df['NES'], errors='coerce').fillna(0)
        direction = np.where(nes > 0, 'up', np.where(nes < 0, 'down', 'neutral'))

        # count: primarily from core_enrichment gene list length; fallback to setSize
        def _count_from_core(x):
            if isinstance(x, str) and x.strip():
                return len([p for p in x.split('/') if p])
            return np.nan

        if 'core_enrichment' in df.columns:
            counts = df['core_enrichment'].map(_count_from_core)
        else:
            counts = pd.Series([np.nan] * len(df))

        if 'setSize' in df.columns:
            counts = counts.fillna(pd.to_numeric(df['setSize'], errors='coerce'))

        counts = counts.fillna(0).astype(int)

        # Override/assign treatment from the CLI-provided label (preserves original intention)
        df_out = pd.DataFrame({
            'treatment' : group,
            'KEGG_ID'   : df['KEGG_ID'],
            'SimpleDesc': simple_desc,
            'direction' : direction.astype(str),
            'count'     : counts
        })

        # Keep only the columns expected downstream
        df_out = df_out[['treatment','KEGG_ID','SimpleDesc','direction','count']]
        all_dfs.append(df_out)

    if not all_dfs:
        sys.exit("Error: No data files were processed.")

    df = pd.concat(all_dfs, ignore_index=True)

    # Group & sum counts for each (treatment, KEGG_ID, SimpleDesc, direction)
    agg = df.groupby(['treatment', 'KEGG_ID', 'SimpleDesc', 'direction'], as_index=False)['count'].sum()

    # Collect the unique treatments in the order they were specified
    treatments_ordered = [t for (_, t) in args.data]
    # Remove duplicates while preserving order
    treatments_unique = list(dict.fromkeys(treatments_ordered))

    # Build a complete grid of (KEGG_ID, treatment) so missing combos become 0
    kegg_ids = agg['KEGG_ID'].unique().tolist()
    all_combos = pd.DataFrame(
        list(itertools.product(kegg_ids, treatments_unique)),
        columns=['KEGG_ID', 'treatment']
    )
    merged_df = pd.merge(all_combos, agg, on=['KEGG_ID','treatment'], how='left')

    # Fill missing
    merged_df['direction'] = merged_df['direction'].fillna("neutral")
    merged_df['count'] = merged_df['count'].fillna(0)
    merged_df['SimpleDesc'] = merged_df['SimpleDesc'].fillna("")

    # We want to sort pathways alphabetically by SimpleDesc alone
    # (No fallback to KEGG_ID).
    # Build a map KEGG_ID -> first SimpleDesc encountered
    desc_map = agg.groupby('KEGG_ID')['SimpleDesc'].first().fillna('').to_dict()

    # Sort KEGG_ID by alphabetical order of its SimpleDesc
    def sort_key(k):
        return desc_map[k].lower()

    sorted_kegg = sorted(kegg_ids, key=sort_key)

    # Convert 'treatment' to categorical so we can sort by the user order
    merged_df['treatment'] = pd.Categorical(
        merged_df['treatment'],
        categories=treatments_unique,
        ordered=True
    )

    # Sort by (KEGG_ID, treatment) in the final DataFrame
    merged_df = merged_df.sort_values(
        by=['KEGG_ID', 'treatment'],
        key=lambda col: col.map(sort_key) if col.name == 'KEGG_ID' else col
    ).reset_index(drop=True)

    return merged_df, treatments_unique, sorted_kegg, desc_map


# =====================
# (1) SINGLE PLOT: Diverging barplot with symmetrical log scale (no +1)
# =====================
def log_diverge(count, direction):
    """
    If count=0 => 0
    Else => ±log10(count)
    """
    if count <= 0:
        return 0
    val = math.log10(count)
    if direction == 'up':
        return val
    elif direction == 'down':
        return -val
    else:
        return 0

def log_diverge_formatter(x, pos):
    """
    Invert symmetrical log scale back to original gene count.
      If x=0 => label "0"
      If x>0 => label = 10^x (rounded)
      If x<0 => label = -(10^abs(x)) (rounded)
    """
    if abs(x) < 1e-12:
        return "0"
    mag = 10**(abs(x))
    mag_int = int(round(mag))
    return f"-{mag_int}" if x < 0 else f"{mag_int}"

def make_single_plot(merged_df, treatments_unique, sorted_kegg, desc_map, output_filename):
    """
    Create a single diverging barplot with symmetrical log scale.
    """
    # Y positions for each KEGG pathway with a bigger gap
    bar_height = 0.15
    gap = 0.1
    y_positions = []
    group_centers = {}
    current_y = 0

    for kegg in sorted_kegg:
        subset = merged_df[merged_df['KEGG_ID'] == kegg]
        idx_list = subset.index.tolist()
        these_ys = []
        for _ in idx_list:
            these_ys.append(current_y)
            current_y += bar_height
        group_centers[kegg] = np.mean(these_ys)
        current_y += gap
        y_positions.extend(these_ys)

    merged_df['y_pos'] = y_positions

    fig, ax = plt.subplots(figsize=(15, max(6, current_y)))

    # Plot bars
    for i, row in merged_df.iterrows():
        tmt = row['treatment']
        up_col, down_col = TREATMENT_COLORS.get(tmt, ("gray", "lightgray"))
        val = log_diverge(row['count'], row['direction'])
        color = "grey"
        if row['direction'] == 'up':
            color = up_col
        elif row['direction'] == 'down':
            color = down_col
        ax.barh(row['y_pos'], val, color=color, height=bar_height, align='edge')

    # Y ticks at group centers
    tick_positions = [group_centers[k] for k in sorted_kegg]
    tick_labels = [desc_map[k] for k in sorted_kegg]
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
    # 🔹 Enlarge KEGG pathway names (y-axis) and x-axis tick labels by 120%
    for lab in ax.get_yticklabels() + ax.get_xticklabels():
        lab.set_fontsize(lab.get_size() * 1.2)


    # X-axis ticks
    all_vals = [log_diverge(r['count'], r['direction']) for _, r in merged_df.iterrows()]
    xmin, xmax = min(all_vals), max(all_vals)
    left_bound = int(math.floor(xmin))
    right_bound = int(math.ceil(xmax))
    x_ticks = range(left_bound, right_bound + 1)

    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(log_diverge_formatter))
    ax.set_xlabel("Number of Genes (Log₁₀ scale)", fontsize=12)

    ax.set_title("Diverging Barplot of KEGG Pathways by Treatment (Log Scale)")

    # Vertical line at x=0
    ax.axvline(0, color='black', linewidth=1)

    # Dashed horizontal lines between pathways
    for kegg in sorted_kegg[1:]:
        y_boundary = merged_df.loc[merged_df['KEGG_ID'] == kegg, 'y_pos'].min() - (gap / 2)
        ax.axhline(y_boundary, color='black', linestyle='dashed', linewidth=0.5, alpha=0.8)

    ax.set_ylim(merged_df['y_pos'].min(), merged_df['y_pos'].max() + bar_height)

    # Single legend with split‐color boxes
    legend_handles = []
    for tmt in treatments_unique:
        up_col, down_col = TREATMENT_COLORS.get(tmt, ("gray", "lightgray"))
        patch = SplitColorPatch(up_col, down_col, label=tmt)
        legend_handles.append(patch)

    ax.legend(
        handles=legend_handles,
        labels=treatments_unique,
        title="Treatment\n(Up=dark, Down=light)",
        handler_map={SplitColorPatch: HandlerSplitColorPatch()},
        loc='upper right'
    )

    plt.tight_layout()
    plt.savefig(output_filename, dpi=600)  # Double resolution
    plt.close(fig)
    print(f"Single diverging barplot saved to {output_filename}")

# =====================
# (2) SUBPLOTS FIGURE: Bubble chart
#  - x=0 => down, x=1 => up
#  - size of bubble => 'count'
#  - sharey=True so pathways align
#  - y-axis labels only on the first subplot
#  - show the left vertical axis line (spine) on ALL subplots
#  - remove the horizontal lines connecting the ticks to the spine in subplots > 0
#  - gray horizontal grid lines in all subplots
#  - legend "Gene Count" on the right side
# =====================
def make_subplots_bubble_plot(merged_df, treatments_unique, sorted_kegg, desc_map, output_filename):
    """
    One subplot per treatment, each showing two columns: x=0 for down, x=1 for up.
    The bubble area is proportional to the gene count.
    We share the same y-axis (pathway order).
    """
    # Build a map from KEGG_ID -> index (alphabetical by SimpleDesc)
    pathway_index = {kegg: i for i, kegg in enumerate(sorted_kegg)}
    num_pathways = len(sorted_kegg)

    # Figure sizing
    fig_width = max(6, (2/3) * len(treatments_unique) * 5)
    fig_height = max(6, num_pathways * 0.3)

    fig, axes = plt.subplots(
        ncols=len(treatments_unique),
        figsize=(fig_width, fig_height),
        sharey=True
    )

    if len(treatments_unique) == 1:
        axes = [axes]

    # Max count for bubble scaling
    max_count = merged_df['count'].max()
    scale_factor = 400 / max_count if max_count > 0 else 1

    for i, (ax, tmt) in enumerate(zip(axes, treatments_unique)):
        sub = merged_df[merged_df['treatment'] == tmt].copy()

        # Horizontal grid lines in every subplot
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='gray', alpha=0.5, linestyle='--', linewidth=0.3)

        # Ensure the left spine is visible on ALL subplots
        ax.spines['left'].set_visible(True)

        # Plot bubbles
        for idx, row in sub.iterrows():
            c = row['count']
            d = row['direction']
            y = pathway_index[row['KEGG_ID']]
            up_col, down_col = TREATMENT_COLORS.get(tmt, ("gray", "lightgray"))

            if c <= 0:
                continue

            if d == 'up':
                ax.scatter(1, y, s=c * scale_factor, color=up_col, alpha=0.8)
            elif d == 'down':
                ax.scatter(0, y, s=c * scale_factor, color=down_col, alpha=0.8)
            # neutral => skip

        ax.set_xlim(-0.5, 1.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Down", "Up"])   # no fontsize here
        ax.set_title(tmt, fontsize=12)

        # For the first subplot:
        if i == 0:
            # Show y-axis tick labels
            ax.set_yticks(range(num_pathways))
            labels = [desc_map[k] for k in sorted_kegg]
            ax.set_yticklabels(labels)       # no fontsize here
        else:
            # Hide y-axis labels for other subplots
            ax.tick_params(labelleft=False)
            ax.tick_params(axis='y', which='both', length=0)

        # 🔹 Enlarge KEGG pathway names (y-axis) and x-axis tick labels by 120%
        for lab in ax.get_yticklabels() + ax.get_xticklabels():
            lab.set_fontsize(lab.get_size() * 1.2)


    # Add a bubble size legend on the right side center
    if max_count > 0:
        # Always show "nice" breakpoints up to and including the first value above max_count
        candidates = [1, 5, 10, 25, 50, 100, 200, 400, 800, 1600, 3200]
        legend_counts = []
        for c in candidates:
            if c <= max_count:
                legend_counts.append(c)
            else:
                legend_counts.append(c)
                break

        # Compute a separate scale so the largest legend bubble is area=400
        legend_max = legend_counts[-1]
        legend_scale = 400.0 / legend_max

        handles = []
        labels = []
        for lc in legend_counts:
            size = lc * legend_scale
            h = axes[-1].scatter([], [], s=size, color='gray', alpha=0.7)
            handles.append(h)
            labels.append(str(lc))

        axes[-1].legend(
            handles,
            labels,
            title="Gene Count",
            loc='center left',
            bbox_to_anchor=(1.05, 0.5),
            labelspacing=1.2       # increase vertical space between entries
        )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_filename, dpi=600)  # Double resolution
    plt.close(fig)
    print(f"Subplot bubble figure saved to {output_filename}")

def _nice_round_numbers_up_to(n):
    """
    Return a list of "nice" round numbers up to n.
    E.g. if n=120, might return [1, 5, 10, 50, 100].
    """
    candidates = [1, 5, 10, 25, 50, 100, 200, 500, 1000, 5000, 10000]
    out = []
    for c in candidates:
        if c <= n:
            out.append(c)
        else:
            break
    return out

def main():
    args = parse_args()

    merged_df, treatments_unique, sorted_kegg, desc_map = read_and_aggregate(args)

    # (1) Single plot with symmetrical log scale
    single_plot_out = args.o
    make_single_plot(merged_df, treatments_unique, sorted_kegg, desc_map, single_plot_out)

    # (2) Subplot bubble chart
    base, ext = os.path.splitext(args.o)
    subplot_out = f"{base}_subplots{ext}"
    make_subplots_bubble_plot(merged_df, treatments_unique, sorted_kegg, desc_map, subplot_out)

if __name__ == "__main__":
    main()
