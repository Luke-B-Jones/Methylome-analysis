#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
make_operon_promoters.py  (RNA-aware)

Compute 600 bp upstream promoter regions for:
  1) Operons given in -i (contig start end strand operon_id), and
  2) Any GFF rows whose COLUMN 2 == 'RNA' from -f (treated like operons).

Promoters are truncated to avoid intruding into:
  - the nearest approach-side interval from the union of:
      * all operon bodies, and
      * all obstructing GFF features (genes, pseudogenes, RNAs, etc.; excluding 'region').

Coordinates:
  - Inputs: 0-based half-open (operons), GFF is 1-based inclusive (converted here).
  - Output: BED6 (0-based, half-open).

Strand rules:
  + strand: promoter = [start-600, start)
  - strand: promoter = [end, end+600)
  Circular wrap handled (split into ≤2 segments).

CLI:
  python make_operon_promoters.py -i operons.tsv -fai genome.fasta.fai -f features.gff -o upstream.bed
"""

import argparse
import csv
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import pandas as pd
import bisect
import re

PROM_LEN = 600  # promoter length in bp


@dataclass
class Target:
    contig: str
    start: int
    end: int
    strand: str
    name: str         # operon_id for operons; Name/locus_tag for RNA
    kind: str         # 'operon' or 'rna'


def load_fai(fai_path: str) -> Dict[str, int]:
    lengths = {}
    with open(fai_path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig = parts[0]
            try:
                L = int(parts[1])
            except ValueError:
                continue
            lengths[contig] = L
    if not lengths:
        raise ValueError(f"No contig lengths parsed from {fai_path}")
    return lengths


def load_operon_targets(tsv_path: str) -> List[Target]:
    """
    Load operons from a 5-column TSV: contig, start, end, strand, operon_id.
    Enforce exactly 5 columns.
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t", header=None, dtype=str, comment="#")
    except Exception:
        df = pd.read_csv(tsv_path, delim_whitespace=True, header=None, dtype=str, comment="#")

    if df.shape[1] != 5:
        raise ValueError(
            f"{tsv_path}: expected exactly 5 columns (contig start end strand operon_id), got {df.shape[1]}"
        )
    df.columns = ["contig", "start", "end", "strand", "operon_id"]

    df["start"] = pd.to_numeric(df["start"], errors="raise").astype(int)
    df["end"]   = pd.to_numeric(df["end"],   errors="raise").astype(int)
    df["strand"] = df["strand"].str.strip().map(lambda s: s[0] if s else "")
    if not df["strand"].isin(["+", "-"]).all():
        bad = df.loc[~df["strand"].isin(["+", "-"])].head(5).to_dict(orient="records")
        raise ValueError(f"Strand column must be '+' or '-'. Bad examples: {bad}")

    targets = [
        Target(
            contig=r.contig,
            start=int(r.start),
            end=int(r.end),
            strand=r.strand,
            name=str(r.operon_id),     # BED name for operons = bare operon_id
            kind="operon",
        )
        for r in df.itertuples(index=False)
    ]
    return targets


def parse_attrs(attrs: str) -> Dict[str, str]:
    d = {}
    for item in attrs.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def load_obstructing_features_gff3(gff_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse GFF3 and return per-contig 0-based half-open intervals for all features
    EXCEPT 'region' rows. (Used to clip promoters.)
    """
    intervals_by_contig: Dict[str, List[Tuple[int, int]]] = {}
    with open(gff_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, source, feature, start, end, *_rest = parts
            if feature.lower() == "region":
                continue
            try:
                start1 = int(start)
                end1   = int(end)
            except ValueError:
                continue
            # GFF: 1-based inclusive -> 0-based half-open
            start0 = max(0, start1 - 1)
            end0   = max(start0, end1)
            intervals_by_contig.setdefault(contig, []).append((start0, end0))
    # sort & merge per contig
    for contig, ivs in intervals_by_contig.items():
        if not ivs:
            continue
        ivs.sort()
        merged: List[Tuple[int, int]] = []
        s0, e0 = ivs[0]
        for s, e in ivs[1:]:
            if s <= e0:
                e0 = max(e0, e)
            else:
                merged.append((s0, e0))
                s0, e0 = s, e
        merged.append((s0, e0))
        intervals_by_contig[contig] = merged
    return intervals_by_contig


def load_rna_targets_from_gff(gff_path: str) -> List[Target]:
    """
    Read GFF3; any row with COLUMN 2 == 'RNA' is treated as a target (like an operon).
    Convert to 0-based half-open. Name preference: Name=, else locus_tag=, else RNA_<start>.
    """
    targets: List[Target] = []
    with open(gff_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, source, feature, start, end, score, strand, phase, attrs = parts
            if source != "RNA":
                continue
            try:
                start1 = int(start)
                end1   = int(end)
            except ValueError:
                continue
            start0 = max(0, start1 - 1)
            end0   = max(start0, end1)
            strand = strand.strip()
            if strand not in {"+", "-"}:
                # If RNA strand missing/invalid, skip (promoters are strand-specific)
                continue
            amap = parse_attrs(attrs)
            name = amap.get("Name") or amap.get("locus_tag") or f"RNA_{start0}"
            targets.append(
                Target(contig=contig, start=start0, end=end0, strand=strand, name=name, kind="rna")
            )
    return targets


def wrap_segments(L: int, raw_start: int, raw_end: int) -> List[Tuple[int, int]]:
    """
    Convert possibly out-of-bounds [raw_start, raw_end) into one or two segments within [0, L).
    """
    segs: List[Tuple[int, int]] = []
    if raw_start < 0:
        segs.append((L + raw_start, L))
        segs.append((0, max(0, raw_end)))
    elif raw_end > L:
        segs.append((max(0, raw_start), L))
        segs.append((0, raw_end - L))
    else:
        segs.append((max(0, raw_start), min(L, raw_end)))
    return [(s, e) for (s, e) in segs if e > s]


def deoverlap_nonwrapping(segments):
    """
    Strand-agnostic, coordinate-space de-overlap.

    Accepts tuples with >=4 fields:
      (start, end, contig, name, *extras)

    - Sort by (contig, start, end).
    - Ensure each segment starts >= previous end on that contig.
    - Preserve any extra trailing fields (e.g., id/kind/strand).
    """
    if not segments:
        return segments

    START, END, CONTIG = 0, 1, 2
    segments = sorted(segments, key=lambda x: (x[CONTIG], x[START], x[END]))
    out = []
    prev_end_by_contig: Dict[str, int] = {}

    for seg in segments:
        s, e, contig = seg[START], seg[END], seg[CONTIG]
        prev_end = prev_end_by_contig.get(contig, -1)

        if s < prev_end:
            s = prev_end

        if e > s:
            seg = (s, e, contig, *seg[3:])
            out.append(seg)
            prev_end_by_contig[contig] = e

    return out


def build_intervals_from_targets(targets: List[Target]) -> Dict[str, List[Tuple[int, int]]]:
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for t in targets:
        by_contig.setdefault(t.contig, []).append((t.start, t.end))
    for contig, ivs in by_contig.items():
        ivs.sort()
    return by_contig


def _sorted_starts_ends(intervals: List[Tuple[int, int]]) -> Tuple[List[int], List[int]]:
    starts = [s for s, _ in intervals]
    ends   = [e for _, e in intervals]
    starts.sort()
    ends.sort()
    return starts, ends


def nearest_limit_for_target(
    t: Target,
    intervals_by_contig: Dict[str, List[Tuple[int, int]]],
) -> Optional[int]:
    """
    For '+' strand, we approach from the LEFT: clamp start by the nearest END <= TSS (= start).
    For '−' strand, we approach from the RIGHT: clamp end by the nearest START >= TSS (= end).

    Return a single integer limit to apply on the approach side, or None if no limit.
    """
    intervals = intervals_by_contig.get(t.contig, [])
    if not intervals:
        return None
    starts, ends = _sorted_starts_ends(intervals)
    if t.strand == "+":
        tss = t.start
        idx = bisect.bisect_right(ends, tss) - 1
        return ends[idx] if idx >= 0 else None
    else:
        tss = t.end
        idx = bisect.bisect_left(starts, tss)
        return starts[idx] if idx < len(starts) else None


def build_promoters(
    targets: List[Target],
    lengths: Dict[str, int],
    features_by_contig: Optional[Dict[str, List[Tuple[int, int]]]] = None,
) -> List[Tuple[str, int, int, str, int, str]]:
    """
    Returns BED6 rows: (contig, start, end, name, score, strand)
    """
    # Bodies of targets (operons + RNAs) as obstructions too
    body_intervals = build_intervals_from_targets(targets)

    # Combine body intervals + all features as obstructions
    combined: Dict[str, List[Tuple[int, int]]] = {k: v[:] for k, v in body_intervals.items()}
    if features_by_contig:
        for contig, ivs in features_by_contig.items():
            combined.setdefault(contig, []).extend(ivs)
    # Merge
    for contig, ivs in combined.items():
        if not ivs:
            continue
        ivs.sort()
        merged: List[Tuple[int, int]] = []
        s0, e0 = ivs[0]
        for s, e in ivs[1:]:
            if s <= e0:
                e0 = max(e0, e)
            else:
                merged.append((s0, e0))
                s0, e0 = s, e
        merged.append((s0, e0))
        combined[contig] = merged

    # Build promoters per target
    pieces_by_strand: Dict[str, List[Tuple[int, int, str, str, str]]] = {"+": [], "-": []}

    for t in targets:
        L = lengths.get(t.contig)
        if L is None:
            raise ValueError(f"Contig {t.contig!r} not found in .fai")

        # Compute raw promoter span
        if t.strand == "+":
            raw_start = t.start - PROM_LEN
            raw_end   = t.start
            limit = nearest_limit_for_target(t, combined)
            if limit is not None:
                raw_start = max(raw_start, limit)
        else:
            raw_start = t.end
            raw_end   = t.end + PROM_LEN
            limit = nearest_limit_for_target(t, combined)
            if limit is not None:
                raw_end = min(raw_end, limit)

        segs = wrap_segments(L, raw_start, raw_end)
        for s, e in segs:
            pieces_by_strand[t.strand].append((s, e, t.contig, t.name, t.kind))

    # De-overlap per strand (strand-agnostic in coordinate space)
    for strand in ["+", "-"]:
        pieces_by_strand[strand] = deoverlap_nonwrapping(pieces_by_strand[strand])

    # Emit BED6
    rows: List[Tuple[str, int, int, str, int, str]] = []
    for strand in ["+", "-"]:
        for s, e, contig, name, _kind in pieces_by_strand[strand]:
            if e > s:
                rows.append((contig, s, e, name, 0, strand))

    # Sort by contig, start for stability (mixed IDs)
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    return rows


def write_bed(rows: List[Tuple[str, int, int, str, int, str]], out_path: str):
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for contig, s, e, name, score, strand in rows:
            w.writerow([contig, s, e, name, score, strand])


def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute 600 bp upstream promoter regions for operons AND any GFF rows with column2=='RNA'."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Operon TSV with EXACT columns: contig  start  end  strand  operon_id (0-based).")
    ap.add_argument("-fai", "--fai", required=True,
                    help="FASTA .fai index giving contig lengths.")
    ap.add_argument("-f", "--features", required=True,
                    help="GFF3 of genome features. Any row with COLUMN 2 == 'RNA' is also used as a target.")
    ap.add_argument("-o", "--out", required=True,
                    help="Output BED6 file of upstream promoters.")
    return ap.parse_args()


def enforce_max_overlap(
    rows: List[Tuple[str, int, int, str, int, str]],
    max_overlap: int = 50,
) -> List[Tuple[str, int, int, str, int, str]]:
    """
    Ensure that ANY overlap between promoters on the same contig is <= max_overlap.
    If two adjacent promoters overlap by > max_overlap, trim BOTH:
      - Find the meeting point (midpoint between the original boundary)
      - Keep 25 bp on each side when max_overlap == 50 (generalized for any even/odd value)
    Operates in coordinate space only (strand-agnostic).

    Input/Output rows are BED6 tuples: (contig, start, end, name, score, strand)
    """

    if not rows:
        return rows

    # group by contig
    by_contig: Dict[str, List[List]] = {}
    for contig, s, e, name, score, strand in rows:
        by_contig.setdefault(contig, []).append([contig, int(s), int(e), name, score, strand])

    half_left = max_overlap // 2
    half_right = max_overlap - half_left  # handles odd max_overlap gracefully

    adjusted: List[Tuple[str, int, int, str, int, str]] = []

    for contig, items in by_contig.items():
        # sort by start then end
        items.sort(key=lambda r: (r[1], r[2]))

        # single left->right pass; trimming current's end and next's start when needed
        for i in range(len(items) - 1):
            curr = items[i]
            nxt  = items[i + 1]
            # current interval: [curr[1], curr[2]), next: [nxt[1], nxt[2])
            overlap = curr[2] - nxt[1]
            if overlap > max_overlap:
                # meeting point (integer midpoint between the two original boundaries)
                meeting = (curr[2] + nxt[1]) // 2
                # new boundaries to keep exactly max_overlap (25/25 when 50)
                new_end_curr  = meeting + half_right
                new_start_nxt = meeting - half_left

                # enforce monotonicity/minimal length
                # don't move beyond original extents (we only shrink)
                new_end_curr  = min(new_end_curr, curr[2])
                new_start_nxt = max(new_start_nxt, nxt[1])

                # ensure positive lengths; let 1 bp minimum survive if heavily trimmed
                if new_end_curr <= curr[1]:
                    new_end_curr = curr[1] + 1
                if nxt[2] <= new_start_nxt:
                    new_start_nxt = nxt[2] - 1
                    if new_start_nxt < nxt[1]:
                        new_start_nxt = nxt[1]  # fallback; may produce 0/1 bp region

                curr[2] = max(curr[1] + 1, new_end_curr)
                nxt[1]  = min(nxt[2] - 1, new_start_nxt)

        # filter out any degenerate intervals (start >= end) after trimming
        for rec in items:
            if rec[2] > rec[1]:
                adjusted.append(tuple(rec))  # type: ignore

    # final stable sort for output
    adjusted.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    return adjusted


def main():
    args = parse_args()
    lengths = load_fai(args.fai)

    # Targets: operons from -i plus RNA entries from -f (column 2 == 'RNA')
    operon_targets = load_operon_targets(args.input)
    rna_targets    = load_rna_targets_from_gff(args.features)
    all_targets    = operon_targets + rna_targets
    if not all_targets:
        raise SystemExit("No targets found (neither operons nor RNA entries).")

    # Obstructions for clipping: union of target bodies + all GFF features (except 'region')
    features_by_contig = load_obstructing_features_gff3(args.features)

    # Build raw promoters (with wrap + truncation against obstructions)
    rows = build_promoters(all_targets, lengths, features_by_contig)

    # Enforce promoter–promoter overlap ≤ 50 bp (25 bp each side of meeting point)
    rows = enforce_max_overlap(rows, max_overlap=50)

    # Write BED6
    write_bed(rows, args.out)

if __name__ == "__main__":
    main()
