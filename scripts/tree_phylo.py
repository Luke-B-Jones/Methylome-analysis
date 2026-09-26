#!/usr/bin/env python3
"""
Peer‑review‑ready pipeline for visualising relationships among DNA/RNA methyl‑transferases
========================================================================================
This version keeps every feature from *tree_fixed.py* and adds reviewer‑level robustness:
  • high‑quality MSA + optional trimAl
  • model‑selected IQ‑TREE with support labels
  • BLOSUM62‑weighted or MSA‑derived %ID
  • seaborn clustered heat‑maps
  • optional SVG export & temp retention

Usage (as before):
python mtase_phylo_heatmap.py -i names.tsv -fa sequences.faa -re refs.fna \
                              -o1 tree1.png -o2 tree2.png --hm1 hm1.png --hm2 hm2.png
Additional switches: --no-trim  --svg  --msa-id  --keep-temp
"""

import argparse
import sys
import os
import tempfile
import subprocess
import shutil
from typing import Dict, List, Sequence

from PIL import Image
# Allow very large scientific figures (trees with many taxa)
Image.MAX_IMAGE_PIXELS = None


from Bio import SeqIO, pairwise2

FONT_SCALE = 1.1  # scale factor for labels and keys

# ----------------------------------------------------------
# Substitution matrix – honour new Biopython location first
# ----------------------------------------------------------
try:
    from Bio.Align import substitution_matrices
    BLOSUM62 = substitution_matrices.load("BLOSUM62")
except Exception:  # pragma: no cover – fall‑back for older versions
    from Bio.SubsMat import MatrixInfo as _MI  # noqa: WPS433 – legacy import
    BLOSUM62 = _MI.blosum62

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Optional heat‑map improvement
try:
    import seaborn as sns  # type: ignore
    SEABORN = True
except ImportError:
    SEABORN = False

# Tree rendering
try:
    from ete3 import Tree, TreeStyle, NodeStyle, faces  # type: ignore
except ImportError:
    sys.exit("Error: ete3 not installed.  pip install ete3")

# Ensure a Qt context for ETE (head‑less safe)
try:
    from PyQt5.QtGui import QGuiApplication  # type: ignore
except ImportError:
    try:
        from PyQt4.QtGui import QGuiApplication  # type: ignore
    except ImportError:
        QGuiApplication = None

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------
def save_colorbar(cmap="viridis", vmin=65, vmax=100,
                  label="Sequence identity (%)",
                  out_png="colorbar.png",
                  orientation="vertical",
                  dpi=1600, fontsize=18, thickness=1.0, length=6.0):
    """
    Save a standalone colorbar as a PNG.
    - thickness: width in inches (if vertical) or height (if horizontal).
    - length: height in inches (if vertical) or width (if horizontal).
    """
    import matplotlib as mpl

    fig, ax = plt.subplots(figsize=(thickness if orientation=="vertical" else length,
                                    length if orientation=="vertical" else thickness))

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=plt.get_cmap(cmap),
                                   norm=norm, orientation=orientation)
    cb.set_label(label, fontsize=fontsize)
    cb.ax.tick_params(labelsize=fontsize*0.9)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved colorbar: {out_png}")


def is_dark_color(hex_color: str) -> bool:
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return (0.299 * r + 0.587 * g + 0.114 * b) < 128


def generate_color_palette(n: int) -> List[str]:
    base = [
        "#e46382",
        "#87d090",
        "#ebd85d",
        "#B2D9EE",
        "#eeb186",
        "#b186be",
        "#84d0d0",
    ]
    return [base[i % len(base)] for i in range(n)]

# -----------------------------------------------------------------------------
# Phylogeny drawing
# -----------------------------------------------------------------------------

def _add_white_border(src: str, dest: str, border: int = 20) -> None:
    """Add a uniform white border around *src* and save to *dest*."""
    img = Image.open(src)
    bg = Image.new("RGB", (img.width + 2 * border, img.height + 2 * border), "white")
    bg.paste(img, (border, border))
    bg.save(dest)

    """Render a cladogram with solid black branches and NO node symbols."""
def style_and_render_tree(
    newick: str,
    labels: Dict[str, str],
    colours: Dict[str, str],
    out_png: str,
    seq_lengths: Dict[str, int],
    italic_labels: Dict[str, bool],
    write_svg: bool = False,
) -> None:


    def custom_layout(node):
        if node.is_leaf():
            label_txt = labels.get(node.name, node.name)
            aa_txt = f"{seq_lengths[node.name]} aa" if node.name in seq_lengths else ""

            # --- Connector line ---
            connector = faces.RectFace(width=20, height=1, fgcolor="black", bgcolor="black")
            connector.opacity = 0.0

            # --- Name face ---
            name_face = faces.TextFace(
                label_txt,
                fsize=23,
                fstyle="italic" if italic_labels.get(node.name, False) else "normal",
            )
            name_face.margin_left = 6
            name_face.margin_right = 6
            if node.name in colours:
                name_face.background.fill = True
                name_face.background.color = colours[node.name]
                name_face.background.opacity = 0.4

            # --- Inline aa face ---
            aa_face = faces.TextFace(aa_txt, fsize=20, fstyle="italic", fgcolor="#333333")
            aa_face.margin_left = 4  # small buffer between name and aa
            aa_face.margin_right = 0

            # --- Add all to branch-right, same column ---
            faces.add_face_to_node(connector, node, column=0, position="branch-right")
            faces.add_face_to_node(name_face, node, column=1, position="branch-right")
            faces.add_face_to_node(aa_face, node, column=2, position="branch-right")

    tree = Tree(newick, format=1)

    # Mid-point root for balance
    mp = tree.get_midpoint_outgroup()
    if mp:
        tree.set_outgroup(mp)

    # Base style: bold black lines
    base = NodeStyle()
    base["shape"] = "none"
    base["size"] = 0
    base["hz_line_width"] = base["vt_line_width"] = 3
    base["hz_line_color"] = base["vt_line_color"] = "black"

    for node in tree.traverse():
        node.set_style(base)

    # Tree style and layout
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.allow_face_overlap = True
    ts.mode = "r"
    ts.scale = 200
    ts.branch_vertical_margin = 12
    ts.show_scale = True
    ts.layout_fn = custom_layout  # Use new layout rendering

    # Render PNG (+ optional SVG)
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False); tmp.close()
    tree.render(tmp.name, w=7000, tree_style=ts, dpi=800)
    _add_white_border(tmp.name, out_png)
    os.remove(tmp.name)
    sys.stderr.write(f"Tree PNG   -> {out_png}\n")

    if write_svg:
        svg = os.path.splitext(out_png)[0] + ".svg"
        tree.render(svg, w=3500, tree_style=ts)
        sys.stderr.write(f"Tree SVG   -> {svg}\n")


# -----------------------------------------------------------------------------
# Heat‑map drawing
# -----------------------------------------------------------------------------
def plot_heatmap(ids: Sequence[str],
                 mat: np.ndarray,
                 label_map: Dict[str, str],
                 italic_map: Dict[str, bool],
                 out_png: str,
                 write_svg: bool = False) -> None:

    labels = [label_map.get(i, i) for i in ids]
    vmin, vmax = 65, 100  # colour scale limits

    plt.figure(figsize=(8, 6))

    if SEABORN:
        ax = sns.heatmap(
            mat,
            cmap="viridis",
            vmin=vmin, vmax=vmax,
            square=True,
            linewidths=0,
            xticklabels=labels,
            yticklabels=labels,
            cbar=False,   # 🔹 disable seaborn's inline colorbar
        )
    else:
        plt.imshow(mat, cmap="viridis", vmin=vmin, vmax=vmax)
        plt.gca().set_aspect("equal")
        # no inline colorbar here either

    plt.xticks(rotation=45, ha="right", fontsize=15*1.2)
    plt.yticks(rotation=0, fontsize=15*1.2)

    for tick, seq_id in zip(plt.gca().get_xticklabels(), ids):
        if italic_map.get(seq_id, False):
            tick.set_fontstyle("italic")

    for tick, seq_id in zip(plt.gca().get_yticklabels(), ids):
        if italic_map.get(seq_id, False):
            tick.set_fontstyle("italic")

    plt.tight_layout()
    plt.savefig(out_png, dpi=800)
    if write_svg:
        plt.savefig(os.path.splitext(out_png)[0] + ".svg")
    plt.close()
    _add_white_border(out_png, out_png)
    sys.stderr.write(f"Heat-map   -> {out_png}\n")

    # 🔹 Save the separate colorbar (bigger, higher res)
    cb_out = os.path.splitext(out_png)[0] + "_colorbar.png"
    save_colorbar(cmap="viridis", vmin=vmin, vmax=vmax,
                  out_png=cb_out, dpi=1600,
                  fontsize=22, thickness=1.5, length=8.0)




# -----------------------------------------------------------------------------
# External tool wrappers
# -----------------------------------------------------------------------------

def _need_prog(name: str) -> None:
    if shutil.which(name) is None:
        sys.exit(f"Error: {name} not found in PATH.")


def check_dependencies() -> None:
    for exe in ("mafft", "iqtree"):
        _need_prog(exe)
    if SEABORN:
        sys.stderr.write("Detected seaborn – clustered heat‑maps enabled.\n")


def run_mafft(in_fasta: str, out_fasta: str) -> None:
    cmd = ["mafft", "--maxiterate", "1000", "--localpair", in_fasta]
    sys.stderr.write("Running MAFFT: " + " ".join(cmd) + "\n")
    res = subprocess.run(cmd, text=True, capture_output=True, check=True)
    with open(out_fasta, "w") as oh:
        oh.write(res.stdout)
    sys.stderr.write(f"Alignment  -> {out_fasta}\n")


def run_trimal(in_aln: str, out_aln: str, skip: bool) -> None:
    if skip or shutil.which("trimal") is None:
        shutil.copy(in_aln, out_aln)
        sys.stderr.write("trimAl not available – using untrimmed alignment.\n")
        return
    cmd = ["trimal", "-automated1", "-in", in_aln, "-out", out_aln]
    sys.stderr.write("Running trimAl: " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    sys.stderr.write(f"Trimmed MSA -> {out_aln}\n")


def run_iqtree(aln: str, pref: str) -> str:
    cmd = ["iqtree", "-s", aln, "-m", "TEST", "-bb", "1000", "-alrt", "1000", "-nt", "AUTO", "-pre", pref]
    sys.stderr.write("Running IQ-TREE: " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    treefile = pref + ".treefile"
    if not os.path.exists(treefile):
        raise RuntimeError(f"IQ-TREE output {treefile} missing.")
    return open(treefile).read().strip()

# -----------------------------------------------------------------------------
# Data wrangling helpers
# -----------------------------------------------------------------------------

def extract_id_fa(rec):
    return rec.id.split("|")[-1] if "|" in rec.id else rec.id


def extract_id_re(rec):
    return rec.id.split("::")[0] if "::" in rec.id else rec.id


def load_records(fasta: str, extractor) -> Dict[str, SeqIO.SeqRecord]:
    d = {}
    for r in SeqIO.parse(fasta, "fasta"):
        k = extractor(r)
        if k not in d:
            r.id = r.name = r.description = k
            d[k] = r
    return d

def parse_label_style(label: str):
    """
    Allow labels wrapped in *...* to be rendered in italics.

    Example:
        *hsdR*  ->  ("hsdR", True)
        hsdR    ->  ("hsdR", False)
    """
    label = label.strip()
    if len(label) >= 2 and label.startswith("*") and label.endswith("*"):
        return label[1:-1], True
    return label, False

def load_names(path: str):
    rows, lbl, italic = [], {}, {}
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue

            parts = [p.strip() for p in ln.split("\t")]
            c1 = parts[0] if len(parts) > 0 else ""
            c2 = parts[1] if len(parts) > 1 else ""
            c3 = parts[2] if len(parts) > 2 else ""

            raw_label = parts[3] if len(parts) > 3 else c1
            label, is_italic = parse_label_style(raw_label)

            assign = parts[-1]
            rows.append((c1, c2, c3, assign))

            if c1:
                lbl[c1] = label
                italic[c1] = is_italic

    return rows, lbl, italic

# -----------------------------------------------------------------------------
# Similarity matrices  —  Biopython-version-agnostic
# -----------------------------------------------------------------------------
def compute_similarity_matrix_pairwise(records: Dict[str, SeqIO.SeqRecord]):
    """
    Return (ids, matrix) where *matrix* contains BLOSUM62-weighted percent
    identities for every sequence pair in *records*.

    The function inspects the runtime signature of
    `pairwise2.align.globaldx` once, selects the right calling convention,
    and caches the resulting helper as `compute_similarity_matrix_pairwise._align`.
    """

    # ----- one-time introspection -------------------------------------------
    if not hasattr(compute_similarity_matrix_pairwise, "_align"):
        import inspect

        sig = inspect.signature(pairwise2.align.globaldx)
        if len(sig.parameters) >= 5:      # old API -> can supply gap penalties
            def _align(a, b, mat, go, ge):
                return pairwise2.align.globaldx(a, b, mat, go, ge)[0]
        else:                             # new API -> gap penalties fixed
            def _align(a, b, mat, *_):
                return pairwise2.align.globaldx(a, b, mat)[0]

        compute_similarity_matrix_pairwise._align = _align  # cache as plain fn

    _align   = compute_similarity_matrix_pairwise._align
    gap_open, gap_extend = -11, -1        # used only for old API

    # ----- build the similarity matrix --------------------------------------
    ids = list(records.keys())
    n   = len(ids)
    mat = np.zeros((n, n), dtype=float)

    for i, id_i in enumerate(ids):
        seq_i = str(records[id_i].seq)
        for j in range(i, n):
            seq_j = str(records[ids[j]].seq)
            aln   = _align(seq_i, seq_j, BLOSUM62, gap_open, gap_extend)

            matches   = sum((a == b) and (a != "-")
                            for a, b in zip(aln.seqA, aln.seqB))
            eff_len   = sum((a != "-") and (b != "-")
                            for a, b in zip(aln.seqA, aln.seqB))
            pid       = 100.0 * matches / eff_len if eff_len else 0.0
            mat[i, j] = mat[j, i] = pid

    return ids, mat


def compute_similarity_matrix_msa(aln: str, ids: List[str]) -> np.ndarray:
    recs = list(SeqIO.parse(aln, "fasta"))
    seq_map = {r.id: str(r.seq) for r in recs}
    aln_len = len(recs[0].seq)
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    for i, a in enumerate(ids):
        for j in range(i, n):
            matches = sum(
                seq_map[a][k] == seq_map[ids[j]][k] and seq_map[a][k] != "-" and seq_map[ids[j]][k] != "-"
                for k in range(aln_len)
            )
            eff = sum(
                seq_map[a][k] != "-" and seq_map[ids[j]][k] != "-" for k in range(aln_len)
            )
            pid = 100.0 * matches / eff if eff else 0.0
            mat[i, j] = mat[j, i] = pid
    return mat

# -----------------------------------------------------------------------------
# CLI parsing
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Build two rectangular cladograms and two heat-maps with optional trimming and SVG export.")

    p.add_argument("-i", "--ids", required=True,
                   help="TSV: col1 fa ID, col2 re ID, col3 e-value, col4 label, last col assignment")
    p.add_argument("-fa", "--fasta", required=True, help="FASTA for col1 IDs (primary).")
    p.add_argument("-re", "--ref", required=True, help="FASTA for col2 IDs (reference).")

    p.add_argument("-o1", "--out1", required=True, help="Cladogram-1 PNG.")
    p.add_argument("-o2", "--out2", required=True, help="Cladogram-2 PNG.")
    p.add_argument("--hm1", "--heatmap1", dest="heatmap1", required=True, help="Heat-map PNG for group 1.")
    p.add_argument("--hm2", "--heatmap2", dest="heatmap2", required=True, help="Heat-map PNG for group 2.")

    # Optional enhancements
    p.add_argument("--no-trim", action="store_true", help="Skip alignment trimming even if trimAl is present.")
    p.add_argument("--svg", action="store_true", help="Also write .svg copies alongside every PNG.")
    p.add_argument("--msa-id", action="store_true",
                   help="Derive percent identities from the final (trimmed) MSA rather than pairwise alignment.")
    p.add_argument("--keep-temp", action="store_true", help="Keep temporary working directory for debugging.")

    return p.parse_args()

# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    check_dependencies()

    # Ensure a Qt application exists for ETE; head-less servers sometimes need this
    if QGuiApplication and QGuiApplication.instance() is None:
        QGuiApplication(sys.argv)

    # Input parsing
    rows, label_map, italic_map = load_names(args.ids)
    recs_fa = load_records(args.fasta, extract_id_fa)
    recs_re = load_records(args.ref, extract_id_re)

    combined = {"1": {}, "2": {}}
    link_map: Dict[str, str] = {}
    cnt = {"1": 0, "2": 0}

    for c1, c2, c3, asgn in rows:
        if asgn not in ("1", "2"):
            continue
        if c1 in recs_fa:
            combined[asgn][c1] = recs_fa[c1]
        if c2 and c2 in recs_re:
            combined[asgn][c2] = recs_re[c2]
        if c2:
            gid = c2  # use reference as group key
            link_map[c1] = gid
            link_map[c2] = gid


    # Process each group separately
    for grp, out_png in (("1", args.out1), ("2", args.out2)):
        if not combined[grp]:
            sys.stderr.write(f"No sequences for cladogram {grp}; skipping.\n")
            continue

        tmpdir = tempfile.TemporaryDirectory() if not args.keep_temp else tempfile.TemporaryDirectory()
        tdir = tmpdir.name

        try:
            # Build combined FASTA
            combined_fa = os.path.join(tdir, "combined.fasta")
            SeqIO.write(combined[grp].values(), combined_fa, "fasta")

            # Align and optionally trim
            aln_fa = os.path.join(tdir, "aligned.fasta")
            run_mafft(combined_fa, aln_fa)

            trimmed_fa = os.path.join(tdir, "aligned.trim.fasta")
            run_trimal(aln_fa, trimmed_fa, skip=args.no_trim)

            # Phylogeny
            newick = run_iqtree(trimmed_fa, os.path.join(tdir, "iq"))

            # Colours for paired IDs
            # Build colour mapping PER GROUP, preserving all pairings
            pair_ids = {}
            for n, gid in link_map.items():
                if n in combined[grp]:
                    pair_ids.setdefault(gid, []).append(n)


            gids = sorted(pair_ids)
            palette = generate_color_palette(len(gids))

            colours = {}
            for i, gid in enumerate(gids):
                for n in pair_ids[gid]:
                    colours[n] = palette[i]

            style_and_render_tree(
                newick,
                label_map,
                colours,
                out_png,
                seq_lengths={k: len(v.seq) for k, v in combined[grp].items()},
                italic_labels=italic_map,
                write_svg=args.svg,
            )


            # Heat-map on *primary* IDs only
            primary_seqs = {k: v for k, v in combined[grp].items() if k in recs_fa}
            if not primary_seqs:
                sys.stderr.write(f"Group {grp}: no primary sequences for heat-map.\n")
                continue

            if args.msa_id:
                ids = list(primary_seqs.keys())
                mat = compute_similarity_matrix_msa(trimmed_fa, ids)
            else:
                ids, mat = compute_similarity_matrix_pairwise(primary_seqs)

            hm_png = args.heatmap1 if grp == "1" else args.heatmap2
            plot_heatmap(ids, mat, label_map, italic_map, hm_png, write_svg=args.svg)
        finally:
            if args.keep_temp:
                sys.stderr.write(f"Temporary files kept in {tdir}\n")
            else:
                tmpdir.cleanup()


if __name__ == "__main__":
    main()