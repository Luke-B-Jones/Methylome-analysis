#!/usr/bin/env python3
"""
DOM_STAN_coverage_grid.py

Coverage plots arranged as a treatment x contig matrix:
    columns = treatment/base (H2O2, VANC, AMPI, ...)
    rows    = contigs/regions

Each panel contains ONLY the matched Treatment_* and Control_* groups for that
one treatment + contig combination.  X is absolute genomic position, not
relative position.  Replicates are summarized as a mean line with a variation
ribbon (SD, SEM, or 95% CI).

Example:
    python DOM_STAN_coverage_grid.py \
      -i HT1/coverage_6mA.txt Treatment_H2O2 \
      -i HT2/coverage_6mA.txt Treatment_H2O2 \
      -i HT3/coverage_6mA.txt Treatment_H2O2 \
      -i HC1/coverage_6mA.txt Control_H2O2 \
      -i HC2/coverage_6mA.txt Control_H2O2 \
      -i HC3/coverage_6mA.txt Control_H2O2 \
      --outdir . --smooth 200 --one-figure
"""

import argparse
import math
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


# -----------------------------------------------------------------------------
# Appearance
# Known bases map to (Treatment colour, Control colour).
TREATMENT_COLORS: Dict[str, Tuple[str, str]] = {
    "H2O2": ("#1b9e77", "#8dcfbb"),
    "VANC": ("#7570b3", "#bab8d9"),
    "AMPI": ("#d95f02", "#efb080"),
    "METH": ("#3395de", "#b3d9f2"),
    "SNAP": ("#b58b00", "#f0dc8a"),
}

BASE_ALIASES = {
    # Accept the H202 spelling in the example bash as H2O2.
    "H202": "H2O2",
    "H2O2": "H2O2",
}

PRETTY_BASE = {
    "H2O2": r"H$_2$O$_2$",
    "VANC": "VANC",
    "AMPI": "AMPI",
    "METH": "METH",
    "SNAP": "SNAP",
}

MEAN_LW = 1.35
CONTROL_LS = (0, (4.0, 2.2))
RIBBON_ALPHA = 0.16


# -----------------------------------------------------------------------------
# I/O

def load_depth(path: str) -> pd.DataFrame:
    """Load a three-column samtools-depth-like file: chrom, pos, depth."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "pos", "depth"],
        usecols=[0, 1, 2],
    )
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df = df.dropna(subset=["chrom", "pos", "depth"]).copy()
    df["pos"] = df["pos"].astype(np.int64)
    df.sort_values(["chrom", "pos"], inplace=True, kind="mergesort")
    return df


def load_bed_optional(path: Optional[str]) -> Optional[pd.DataFrame]:
    """Load BED chrom/start/end/[name]."""
    if not path:
        return None
    if not os.path.isfile(path):
        raise SystemExit(f"--bed file not found: {path}")

    with open(path, "r", encoding="utf-8") as fh:
        first = next((line.strip() for line in fh if line.strip()), "")

    has_name = len(first.split()) >= 4 if first else False
    names = ["chrom", "start", "end"] + (["name"] if has_name else [])
    usecols = [0, 1, 2] + ([3] if has_name else [])
    bed = pd.read_csv(path, sep="\t", header=None, names=names, usecols=usecols)
    bed["chrom"] = bed["chrom"].astype(str)
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
    bed = bed.dropna(subset=["chrom", "start", "end"]).copy()
    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["end"].astype(int)
    if "name" not in bed.columns:
        bed["name"] = bed["chrom"]
    else:
        bed["name"] = bed["name"].astype(str)
    return bed[["chrom", "start", "end", "name"]]


def natural_key(text: str):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text)]


def infer_regions_from_depths(
    grouped: Dict[Tuple[str, str], List[pd.DataFrame]]
) -> Optional[pd.DataFrame]:
    """Infer one full-span region per contig across every input."""
    spans: Dict[str, List[int]] = {}
    for dfs in grouped.values():
        for df in dfs:
            for chrom, g in df.groupby("chrom", sort=False):
                lo = int(g["pos"].min())
                hi = int(g["pos"].max())
                if chrom not in spans:
                    spans[chrom] = [lo, hi]
                else:
                    spans[chrom][0] = min(spans[chrom][0], lo)
                    spans[chrom][1] = max(spans[chrom][1], hi)

    if not spans:
        return None

    rows = [
        {"chrom": chrom, "start": lo, "end": hi, "name": chrom}
        for chrom, (lo, hi) in spans.items()
    ]
    rows.sort(key=lambda x: natural_key(x["chrom"]))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])


# -----------------------------------------------------------------------------
# Group parsing / colours

def canonical_base(base: str) -> str:
    key = base.strip().upper().replace("-", "")
    return BASE_ALIASES.get(key, key)


def parse_group(group: str) -> Tuple[str, str]:
    """Return (role, base), where role is treatment/control/other."""
    g = group.strip()
    if "_" not in g:
        return "other", canonical_base(g)

    head, tail = g.split("_", 1)
    head = head.lower()
    if head.startswith("treat"):
        role = "treatment"
    elif head.startswith("control"):
        role = "control"
    else:
        role = "other"
    return role, canonical_base(tail)


def lighten(color: str, amount: float = 0.62) -> str:
    """Blend a colour towards white; amount=0 keeps colour, 1 gives white."""
    rgb = np.asarray(to_rgb(color))
    out = rgb + (1.0 - rgb) * amount
    return to_hex(out)


def colors_for_base(base: str) -> Tuple[str, str]:
    if base in TREATMENT_COLORS:
        return TREATMENT_COLORS[base]

    # Deterministic fallback for an unknown treatment name.
    cmap = plt.get_cmap("tab10")
    idx = sum(ord(ch) for ch in base) % 10
    treatment = to_hex(cmap(idx))
    control = lighten(treatment)
    return treatment, control


def pretty_base(base: str) -> str:
    return PRETTY_BASE.get(base, base)


# -----------------------------------------------------------------------------
# Summaries

def slice_region(df: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    mask = (df["chrom"] == chrom) & (df["pos"] >= start) & (df["pos"] <= end)
    return df.loc[mask, ["pos", "depth"]]


def build_group_matrix(series_list: List[pd.Series]) -> pd.DataFrame:
    """Align replicate depth vectors by genomic position."""
    if not series_list:
        return pd.DataFrame()
    mat = pd.concat(series_list, axis=1)
    mat.columns = [f"rep{i + 1}" for i in range(mat.shape[1])]
    return mat.fillna(0.0).sort_index()


def smooth_series(s: pd.Series, window: int) -> pd.Series:
    if window and window > 1:
        return s.rolling(window=window, min_periods=1, center=True).mean()
    return s


def compute_mean_and_ribbon(
    mat: pd.DataFrame,
    ribbon: str,
    smooth: int,
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    if mat.empty:
        empty = pd.Series(dtype=float)
        return empty, empty, empty

    n = mat.shape[1]
    mean = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=1) if n > 1 else pd.Series(0.0, index=mat.index)

    if ribbon == "sd":
        half = sd
    elif ribbon == "sem":
        half = sd / np.sqrt(max(1, n))
    elif ribbon == "ci95":
        half = 1.96 * sd / np.sqrt(max(1, n))
    else:
        raise SystemExit(f"Unknown --ribbon value: {ribbon}")

    mean_sm = smooth_series(mean, smooth)
    half_sm = smooth_series(half, smooth)
    lower = (mean_sm - half_sm).clip(lower=0.0)
    upper = mean_sm + half_sm
    return mean_sm, lower, upper


def summarize_group_region(
    dfs: List[pd.DataFrame],
    chrom: str,
    start: int,
    end: int,
    ribbon: str,
    smooth: int,
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    reps: List[pd.Series] = []
    for df in dfs:
        sub = slice_region(df, chrom, start, end)
        if sub.empty:
            continue
        reps.append(sub.set_index("pos")["depth"].astype(float))
    return compute_mean_and_ribbon(build_group_matrix(reps), ribbon, smooth)


def global_depth_limits(
    grouped: Dict[Tuple[str, str], List[pd.DataFrame]],
    ymax: Optional[float],
) -> Tuple[float, float]:
    if ymax is not None:
        if ymax <= 0:
            raise SystemExit("--ymax must be > 0")
        return 0.0, float(ymax)

    dmax = 0.0
    for dfs in grouped.values():
        for df in dfs:
            if not df.empty:
                dmax = max(dmax, float(df["depth"].max()))
    return 0.0, max(1.0, dmax * 1.05)


# -----------------------------------------------------------------------------
# Plot helpers

def panel_letter(index: int) -> str:
    """0 -> A, 25 -> Z, 26 -> AA."""
    out = ""
    n = index + 1
    while n:
        n, rem = divmod(n - 1, 26)
        out = chr(65 + rem) + out
    return out


def x_scale_and_label(unit: str) -> Tuple[float, str]:
    if unit == "bp":
        return 1.0, "Genomic position (bp)"
    if unit == "kb":
        return 1e3, "Genomic position (kb)"
    return 1e6, "Genomic position (Mb)"


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="both", labelsize=8.5, width=0.8, length=3.5)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune=None))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune=None))
    ax.grid(axis="y", color="#e8e8e8", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def add_figure_legend(fig, ribbon: str, position: str = "bottom"):
    ribbon_label = {
        "sd": "variation ribbon: ±SD",
        "sem": "variation ribbon: ±SEM",
        "ci95": "variation ribbon: 95% CI",
    }[ribbon]

    handles = [
        Line2D([0], [0], color="#222222", lw=1.6, linestyle="-", label="Treatment mean"),
        Line2D([0], [0], color="#777777", lw=1.4, linestyle=CONTROL_LS, label="Control mean"),
        Patch(facecolor="#777777", alpha=RIBBON_ALPHA, edgecolor="none", label=ribbon_label),
    ]

    if position == "bottom":
        loc = "lower center"
        anchor = (0.5, 0.012)
    else:
        loc = "upper center"
        anchor = (0.5, 0.992)

    leg = fig.legend(
        handles=handles,
        loc=loc,
        bbox_to_anchor=anchor,
        frameon=True,
        fancybox=True,
        framealpha=1.0,
        borderpad=0.65,
        handlelength=2.8,
        columnspacing=1.8,
        fontsize=9,
        ncol=3,
    )
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_edgecolor("#bdbdbd")


def plot_grid(
    regions: pd.DataFrame,
    grouped: Dict[Tuple[str, str], List[pd.DataFrame]],
    bases: List[str],
    outdir: str,
    smooth: int,
    ribbon: str,
    ylim: Tuple[float, float],
    x_unit: str,
    dpi: int,
    fname_prefix: str,
):
    """Rows = contigs/regions; columns = treatments."""
    nrows = len(regions)
    ncols = len(bases)
    if nrows == 0 or ncols == 0:
        raise SystemExit("Nothing to plot: no regions or no Treatment_/Control_ groups found.")

    fig_w = max(8.5, 4.25 * ncols)
    fig_h = max(4.0, 2.35 * nrows + 1.05)

    # sharex='row' means all treatments for the same contig use exactly the same genomic x-axis.
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
        sharex="row",
        sharey=True,
    )

    scale, xlabel = x_scale_and_label(x_unit)
    panel_i = 0

    for r, (_, region) in enumerate(regions.iterrows()):
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])
        row_name = str(region["name"])

        for c, base in enumerate(bases):
            ax = axes[r, c]
            style_axis(ax)
            treatment_color, control_color = colors_for_base(base)

            # Treatment and matched control only: exactly one treatment-contig combo per panel.
            for role, color, linestyle, zorder in (
                ("control", control_color, CONTROL_LS, 2),
                ("treatment", treatment_color, "-", 3),
            ):
                dfs = grouped.get((role, base), [])
                mean, lower, upper = summarize_group_region(
                    dfs=dfs,
                    chrom=chrom,
                    start=start,
                    end=end,
                    ribbon=ribbon,
                    smooth=smooth,
                )
                if mean.empty:
                    continue

                x = mean.index.to_numpy(dtype=float) / scale
                ax.fill_between(
                    x,
                    lower.to_numpy(),
                    upper.to_numpy(),
                    color=color,
                    alpha=RIBBON_ALPHA,
                    linewidth=0,
                    zorder=zorder,
                )
                ax.plot(
                    x,
                    mean.to_numpy(),
                    color=color,
                    lw=MEAN_LW,
                    linestyle=linestyle,
                    zorder=zorder + 1,
                )

            ax.set_xlim(start / scale, end / scale)
            ax.set_ylim(*ylim)

            # Small panel letter, like a publication multi-panel figure.
            ax.text(
                0.018,
                0.94,
                panel_letter(panel_i),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                fontweight="bold",
            )
            panel_i += 1

            # Column headers appear once.
            if r == 0:
                ax.set_title(pretty_base(base), fontsize=15, fontweight="bold", pad=10)

            # A row label appears only once, outside the first-column axis.
            if c == 0:
                ax.annotate(
                    row_name,
                    xy=(0, 0.5),
                    xycoords="axes fraction",
                    xytext=(-0.19, 0.5),
                    textcoords="axes fraction",
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=9.5,
                    fontweight="bold",
                    annotation_clip=False,
                )
                ax.set_ylabel("Depth", fontsize=9.5, labelpad=7)
            else:
                ax.tick_params(labelleft=False)

            # Keep absolute genomic tick labels on every row; only the shared axis title is global.
            ax.tick_params(labelbottom=True)

    # Put the shared x-axis title above the key, with the key beneath the full grid.
    fig.supxlabel(xlabel, fontsize=11, y=0.10)
    add_figure_legend(fig, ribbon, position="bottom")

    # Reserve a clean band below the axes for the shared x-label and legend/key.
    fig.subplots_adjust(
        left=0.105,
        right=0.985,
        bottom=0.20,
        top=0.91,
        wspace=0.08,
        hspace=0.23,
    )

    png = os.path.join(outdir, f"{fname_prefix}.png")
    pdf = os.path.join(outdir, f"{fname_prefix}.pdf")
    fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[info] Saved: {png}")
    print(f"[info] Saved: {pdf}")


def plot_separate_panels(
    regions: pd.DataFrame,
    grouped: Dict[Tuple[str, str], List[pd.DataFrame]],
    bases: List[str],
    outdir: str,
    smooth: int,
    ribbon: str,
    ylim: Tuple[float, float],
    x_unit: str,
    dpi: int,
):
    """Optional non-grid mode: one image for each treatment x contig combination."""
    scale, xlabel = x_scale_and_label(x_unit)

    for _, region in regions.iterrows():
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])
        row_name = str(region["name"])

        for base in bases:
            fig, ax = plt.subplots(figsize=(8.0, 3.8))
            style_axis(ax)
            tcol, ccol = colors_for_base(base)

            for role, color, linestyle, zorder in (
                ("control", ccol, CONTROL_LS, 2),
                ("treatment", tcol, "-", 3),
            ):
                mean, lower, upper = summarize_group_region(
                    grouped.get((role, base), []),
                    chrom,
                    start,
                    end,
                    ribbon,
                    smooth,
                )
                if mean.empty:
                    continue
                x = mean.index.to_numpy(dtype=float) / scale
                ax.fill_between(x, lower, upper, color=color, alpha=RIBBON_ALPHA, linewidth=0)
                ax.plot(x, mean, color=color, lw=MEAN_LW, linestyle=linestyle, zorder=zorder)

            ax.set_xlim(start / scale, end / scale)
            ax.set_ylim(*ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Depth")
            ax.set_title(f"{pretty_base(base)} — {row_name}", fontweight="bold")
            add_figure_legend(fig, ribbon, position="bottom")
            fig.subplots_adjust(left=0.11, right=0.98, bottom=0.24, top=0.88)

            safe_contig = re.sub(r"[^A-Za-z0-9_.-]+", "_", row_name)
            safe_base = re.sub(r"[^A-Za-z0-9_.-]+", "_", base)
            path = os.path.join(outdir, f"{safe_base}__{safe_contig}_coverage.png")
            fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"[info] Saved: {path}")


# -----------------------------------------------------------------------------
# CLI

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Plot mean coverage with a variation ribbon. In --one-figure mode, "
            "rows are contigs and columns are treatment bases; each panel contains "
            "the matched Treatment_* and Control_* groups only."
        )
    )
    ap.add_argument(
        "-i",
        "--input",
        action="append",
        nargs=2,
        metavar=("DEPTH_TXT", "GROUP"),
        required=True,
        help="Repeatable depth file + group, e.g. -i sample.txt Treatment_H2O2",
    )
    ap.add_argument(
        "--bed",
        help="Optional BED: chrom start end [name]. If absent, infer one full span per contig.",
    )
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--smooth", type=int, default=50, help="Centered rolling window; 0/1 disables")
    ap.add_argument(
        "--ribbon",
        choices=["sd", "sem", "ci95"],
        default="sd",
        help="Variation ribbon around the replicate mean (default: sd)",
    )
    ap.add_argument(
        "--one-figure",
        action="store_true",
        help="Make the treatment x contig grid (recommended).",
    )
    ap.add_argument(
        "--x-unit",
        choices=["bp", "kb", "Mb"],
        default="Mb",
        help="Display absolute genomic position in these units (default: Mb).",
    )
    ap.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional fixed y-axis maximum. Default uses the global observed max.",
    )
    ap.add_argument("--dpi", type=int, default=400, help="PNG DPI (default: 400)")
    ap.add_argument(
        "--output-name",
        default="coverage_by_treatment_and_contig",
        help="Base filename for --one-figure output",
    )

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # (role, base) -> replicate dataframes. Dict insertion order preserves CLI order.
    grouped: Dict[Tuple[str, str], List[pd.DataFrame]] = defaultdict(list)
    base_order: List[str] = []

    for path, group_name in args.input:
        if not os.path.isfile(path):
            raise SystemExit(f"Depth file not found: {path}")

        role, base = parse_group(group_name)
        if role == "other":
            raise SystemExit(
                f"Group '{group_name}' is not Treatment_* or Control_*. "
                "Use names such as Treatment_H2O2 / Control_H2O2."
            )

        if base not in base_order:
            base_order.append(base)

        df = load_depth(path)
        if df.empty:
            print(f"[warn] No usable rows in: {path}")
        grouped[(role, base)].append(df)

    if not grouped:
        raise SystemExit("No valid inputs loaded.")

    regions = load_bed_optional(args.bed)
    if regions is None or regions.empty:
        regions = infer_regions_from_depths(grouped)
    if regions is None or regions.empty:
        raise SystemExit("No contigs/regions available to plot.")

    # Helpful diagnostics for an incomplete treatment-control pair.
    for base in base_order:
        nt = len(grouped.get(("treatment", base), []))
        nc = len(grouped.get(("control", base), []))
        print(f"[info] {base}: {nt} treatment replicate(s), {nc} control replicate(s)")
        if nt == 0 or nc == 0:
            print(f"[warn] {base} does not have both Treatment_* and Control_* inputs.")

    ylim = global_depth_limits(grouped, args.ymax)
    print(f"[info] Shared y-axis: {ylim[0]:.3g} .. {ylim[1]:.3g}")
    print(f"[info] Grid: {len(regions)} row(s) x {len(base_order)} column(s)")

    if args.one_figure:
        plot_grid(
            regions=regions,
            grouped=grouped,
            bases=base_order,
            outdir=args.outdir,
            smooth=args.smooth,
            ribbon=args.ribbon,
            ylim=ylim,
            x_unit=args.x_unit,
            dpi=args.dpi,
            fname_prefix=args.output_name,
        )
    else:
        plot_separate_panels(
            regions=regions,
            grouped=grouped,
            bases=base_order,
            outdir=args.outdir,
            smooth=args.smooth,
            ribbon=args.ribbon,
            ylim=ylim,
            x_unit=args.x_unit,
            dpi=args.dpi,
        )


if __name__ == "__main__":
    main()
