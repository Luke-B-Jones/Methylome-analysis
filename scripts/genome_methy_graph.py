#!/usr/bin/env python3
"""
genome_methy_graph.py

Per-contig IGV-style stacked barplots showing total modified bases per bin.
All treatments appear in a single row per contig.
Legend is saved as a separate PNG.

"""

from __future__ import annotations

import argparse
import math
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------- IO --------------------------------------------

def read_fai(path: str) -> Dict[str, int]:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["contig", "length", "o", "lb", "lw"],
        usecols=[0, 1],
    )
    if df.empty:
        raise ValueError("FAI file is empty.")
    return dict(zip(df["contig"], df["length"].astype(int)))


def read_sites(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = {"chrom", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {', '.join(sorted(missing))}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end"]).copy()

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Support single-base or interval features
    df["n_bases"] = (df["end"] - df["start"]).clip(lower=1)

    return df


# ----------------------------- Binning ---------------------------------------

def bin_sites_total_bases(
    sites: pd.DataFrame,
    contig_lengths: Dict[str, int],
    bin_size: int,
) -> pd.DataFrame:
    """
    Bin sites and sum total bases per bin.
    """
    sites = sites[sites["chrom"].isin(contig_lengths)].copy()
    sites["bin"] = (sites["start"] // bin_size).astype(int)

    rows = []
    for chrom, L in contig_lengths.items():
        n_bins = int(math.ceil(L / bin_size))
        sub = sites[sites["chrom"] == chrom]
        if sub.empty:
            continue

        g = sub.groupby("bin")["n_bases"].sum()
        for b, v in g.items():
            if b < n_bins:
                rows.append((chrom, b, int(v)))

    return pd.DataFrame(
        rows, columns=["chrom", "bin", "total_bases"]
    )


# ----------------------------- Plotting --------------------------------------

def plot_per_contig_barplots(
    out_dir: str,
    contig_lengths: Dict[str, int],
    treatments: List[str],
    colors: Dict[str, str],
    binned_by_treatment: Dict[str, pd.DataFrame],
    bin_size: int,
):
    """
    One stacked barplot per contig (NO legend).
    """

    contigs = list(contig_lengths.keys())

    # ----- global y-scale -----
    global_max = 0
    for t in treatments:
        if not binned_by_treatment[t].empty:
            global_max = max(global_max, binned_by_treatment[t]["total_bases"].max())
    global_max = max(1, global_max)

    for contig in contigs:
        L = contig_lengths[contig]
        n_bins = int(math.ceil(L / bin_size))

        x = np.arange(n_bins)
        bottom = np.zeros(n_bins)

        fig_width = max(8, min(16, L / 5000))
        fig, ax = plt.subplots(figsize=(fig_width, 3.5))

        # ---- plot bars ----
        for t in treatments:
            df = binned_by_treatment[t]
            sub = df[df["chrom"] == contig]

            y = np.zeros(n_bins)
            for _, r in sub.iterrows():
                if r["bin"] < n_bins:
                    y[int(r["bin"])] = r["total_bases"]

            ax.bar(
                x,
                y,
                bottom=bottom,
                width=1.0,
                color=colors[t],
                edgecolor="none",
            )
            bottom += y

        # ✅ THIS IS WHERE IT GOES
        ax.set_xlim(0, n_bins - 1)
        ax.margins(x=0)

        # ---- axes & labels ----
        xticks = np.linspace(0, n_bins - 1, num=5, dtype=int)
        ax.set_xticks(xticks)
        ax.set_xticklabels(
            [f"{(b * bin_size) / 1000:.1f}" for b in xticks],
            fontsize=11,
        )

        ax.set_xlabel("Position (kb)", fontsize=13)
        ax.set_ylabel("Total bases (≥20% ΔPL) $\\mathit{per}$ bin", fontsize=13)
        ax.set_ylim(0, global_max * 1.1)

        ax.set_title(
            f"{contig}  (bin = {bin_size} bp, length = {L:,} bp)",
            fontsize=14,
        )

        ax.tick_params(axis="both", labelsize=11)

        fig.tight_layout()
        fig.savefig(
            os.path.join(out_dir, f"{contig}_barplot.png"),
            dpi=600,
        )
        plt.close(fig)



def save_legend_png(
    out_png: str,
    treatments: List[str],
    colors: Dict[str, str],
):
    """
    Save a standalone legend PNG.
    """
    fig, ax = plt.subplots(figsize=(2 + 1.2 * len(treatments), 1.2))
    ax.axis("off")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[t])
        for t in treatments
    ]

    ax.legend(
        handles,
        treatments,
        loc="center",
        frameon=False,
        ncol=len(treatments),
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)


# ----------------------------- CLI -------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Per-contig stacked barplots of methylation signal"
    )
    p.add_argument("--fai", required=True, help="FASTA index (.fai)")
    p.add_argument(
        "-i",
        nargs=3,
        action="append",
        metavar=("FILE", "NAME", "COLOR"),
        required=True,
        help="Input file, treatment name, and color",
    )
    p.add_argument("-b", type=int, required=True, help="Bin size (bp)")
    p.add_argument("-o", required=True, help="Output directory")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.b <= 0:
        raise ValueError("Bin size must be > 0")

    os.makedirs(args.o, exist_ok=True)

    contig_lengths = read_fai(args.fai)

    treatments: List[str] = []
    colors: Dict[str, str] = {}
    binned_by_treatment: Dict[str, pd.DataFrame] = {}

    for path, name, color in args.i:
        if name in treatments:
            raise ValueError(f"Duplicate treatment name: {name}")
        treatments.append(name)
        colors[name] = color

        sites = read_sites(path)
        binned = bin_sites_total_bases(
            sites, contig_lengths, args.b
        )
        binned_by_treatment[name] = binned

    plot_per_contig_barplots(
        out_dir=args.o,
        contig_lengths=contig_lengths,
        treatments=treatments,
        colors=colors,
        binned_by_treatment=binned_by_treatment,
        bin_size=args.b,
    )

    save_legend_png(
        out_png=os.path.join(args.o, "legend.png"),
        treatments=treatments,
        colors=colors,
    )

    print(f"[OK] Per-contig barplots + legend written to {args.o}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
