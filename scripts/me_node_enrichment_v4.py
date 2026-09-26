#!/usr/bin/env python3
"""
me_node_enrichment.py — version 5.0.0

A general, window-level analysis of whether pre-treatment methylation entropy
(ME) is associated with later treatment-responsive methylation change.

Version 5 separates four spatial/statistical questions:

1. SITE-OVERLAPPING ME
   ME from true Modkit entropy windows physically overlapping the exact base.

2. LOCAL +/-500-BP CONTEXT
   ME from non-overlapping windows within +/-500 bp of that base.

3. CALLABLE GENOME BACKGROUND
   ME at eligible callable sites sampled across the same
   modification/contig/strand risk set.

4. LATER RESPONSE
   Either continuous |treatment PL - control PL| or the binary event that a
   site crossed the supplied treatment-response threshold.

Every inclusive analysis is run twice:
- adjusted for baseline PL, depth, ME support, sequence composition,
  overlapping-window count, callability and treatment where relevant;
- unadjusted, retaining only sampling weights and genomic-block uncertainty.

The historical strict case-specific matching, treatment/direction analyses and
sensitivity analyses are retained. Additional binary, reverse-orientation and
direct peak tests retain every threshold-positive site with valid ME.

All formal ME endpoints use ``windows.bedgraph``. ``regions.bed`` is QC-only.
The command-line interface remains compatible with Versions 3/4; each ``--me``
path is a directory containing:

    windows.bedgraph
    regions.bed

Example
-------
python me_node_enrichment_v5.py \
  --metadata NEW_META_report.tsv \
  --fasta genome.fasta \
  --fai genome.fasta.fai \
  --gff gene_annot.gff \
  --operons operon_list.tsv \
  --dmr-root DMR_ROOT \
  --me 4mC CONTROL 4mC_ME/ \
  --me 5mC CONTROL 5mC_ME/ \
  --me 6mA CONTROL 6mA_ME/ \
  -o ME_enrichment_v5

Scientific scope
----------------
With pooled control-ME inputs, this is a genomic-location association analysis.
It can test whether baseline control ME marks sites that later change. It
cannot by itself establish biological-replicate reproducibility,
treatment-induced change in ME, stable epialleles, or coherent whole-genome
methylotypes.

Dependencies
------------
Python >= 3.7, numpy, pandas, matplotlib; SciPy is optional.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import gzip
import hashlib
import json
import logging
import math
import os
import random
import re
import shutil
import sqlite3
import sys
import time
import zipfile
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


VERSION = "5.0.0"
LOG = logging.getLogger("me_node_enrichment_v5")

MISSING = {"", ".", "na", "nan", "n/a", "null", "none"}
DNA_BASES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

MOD_COLUMN_RE = re.compile(r"^(?P<treatment>.+)_(?P<mod>[^_]+)$")
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
SITE_RECORD_RE = re.compile(
    rf"(?P<chrom>[^:;|]+):"
    rf"(?P<start>\d+)-(?P<end>\d+)"
    rf"\((?P<strand>[+.\-])\):"
    rf"(?P<reported_delta>{FLOAT})"
    rf"(?:\[(?P<control_pl>{FLOAT}),(?P<treatment_pl>{FLOAT})\])?"
)
TOKEN_SPLIT_RE = re.compile(r"\s*(?:\|\|\||;)\s*")


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass(frozen=True, order=True)
class SiteKey:
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def position(self) -> int:
        return self.start


@dataclass
class CaseSite:
    mod_type: str
    key: SiteKey
    treatments: set[str] = field(default_factory=set)
    reported_deltas: dict[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    control_pl: dict[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    treatment_pl: dict[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    source_regions: set[str] = field(default_factory=set)

    def delta_values(self) -> list[float]:
        values: list[float] = []
        for treatment in self.treatments:
            controls = self.control_pl.get(treatment, [])
            treated = self.treatment_pl.get(treatment, [])
            for c, t in zip(controls, treated):
                values.append(t - c)
        return values

    @property
    def direction(self) -> str:
        values = [value for value in self.delta_values() if value != 0]
        if not values:
            return "unknown"
        if all(value > 0 for value in values):
            return "gain"
        if all(value < 0 for value in values):
            return "loss"
        return "mixed"

    @property
    def max_abs_delta(self) -> float:
        values = self.delta_values()
        if values:
            return max(abs(value) for value in values)
        reported = [
            abs(value)
            for treatment_values in self.reported_deltas.values()
            for value in treatment_values
        ]
        return max(reported, default=float("nan"))


@dataclass
class DMRFile:
    mod_type: str
    treatment: str
    path: Path
    bit: int
    treatment_index: int


@dataclass
class SiteRecord:
    key: SiteKey
    name: str
    mask: int
    pls: tuple[Optional[float], ...]
    depths: tuple[Optional[float], ...]
    post_pls: tuple[Optional[float], ...]
    post_depths: tuple[Optional[float], ...]
    sampling_weight: float = 1.0
    context: str = "unannotated"


@dataclass
class WindowTrack:
    sample: str
    mod_type: str
    path: Path
    stranded: bool
    chrom: str
    strand: str
    starts: np.ndarray
    ends: np.ndarray
    midpoints: np.ndarray
    entropy: np.ndarray
    support: np.ndarray
    median_width: float
    median_step: float
    max_width: int


@dataclass
class SiteScore:
    mod_type: str
    key: SiteKey
    strand_mode: str
    site_me: float
    site_q75: float
    site_support: float
    overlap_windows: int
    domain_250: float
    domain_500: float
    domain_1000: float
    context_500_excluding_site: float
    local_prominence_500: float
    local_percentile_500: float
    local_context_windows_500: int
    local_window_median: float
    local_window_prominence: float
    control_samples: int

    def endpoint(self, name: str) -> float:
        return float(getattr(self, name))


@dataclass
class CaseAnalysis:
    case: CaseSite
    case_record: SiteRecord
    case_score: SiteScore
    local_control_values: np.ndarray
    local_control_weights: np.ndarray
    local_control_keys: tuple[SiteKey, ...]
    global_control_values: np.ndarray
    global_control_weights: np.ndarray
    global_control_keys: tuple[SiteKey, ...]
    local_effect: float
    local_percentile: float
    global_effect: float
    global_percentile: float
    local_ess: float
    global_ess: float
    local_controls: int
    global_controls: int
    local_within_case_max_standardized_mismatch: float
    global_within_case_max_standardized_mismatch: float
    case_control_pl: float
    case_control_depth: float
    case_me_support: float
    case_gc_fraction: float
    case_target_base_fraction: float
    local_expected_control_pl: float
    local_expected_control_depth: float
    local_expected_me_support: float
    local_expected_gc_fraction: float
    local_expected_target_base_fraction: float
    global_expected_control_pl: float
    global_expected_control_depth: float
    global_expected_me_support: float
    global_expected_gc_fraction: float
    global_expected_target_base_fraction: float


# =============================================================================
# CLI AND UTILITIES
# =============================================================================

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Window-level control-ME localization and regional-enrichment "
            "analysis for treatment-discovered methylation-change sites."
        ),
    )
    parser.add_argument("--metadata", "-m", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--operons", "-u", required=True, type=Path)
    parser.add_argument("--dmr-root", required=True, type=Path)
    parser.add_argument(
        "--me",
        nargs=3,
        action="append",
        required=True,
        metavar=("MOD_TYPE", "CONTROL_LABEL", "ME_DIRECTORY"),
        help=(
            "Repeat for each modification/control input. ME_DIRECTORY must "
            "contain windows.bedgraph and regions.bed."
        ),
    )
    parser.add_argument("-o", "--output-dir", required=True, type=Path)

    parser.add_argument(
        "--min-window-support",
        type=int,
        default=10,
        help="Primary minimum reads per entropy window.",
    )
    parser.add_argument(
        "--min-finite-me-coverage",
        type=float,
        default=0.80,
        help=(
            "Finite-entropy coverage QC warning threshold. Coverage below "
            "this value is reported but does not terminate the analysis; "
            "formal tests require valid ME at retained cases and controls."
        ),
    )
    parser.add_argument(
        "--support-sensitivity",
        nargs="*",
        type=int,
        default=[5, 20],
        help="Additional minimum-support thresholds.",
    )
    parser.add_argument(
        "--local-radius",
        type=int,
        default=2000,
        help="Neighbourhood radius for the local-site test.",
    )
    parser.add_argument(
        "--domain-radius",
        type=int,
        default=500,
        choices=[250, 500, 1000],
        help="Primary radius for the broad-domain test.",
    )
    parser.add_argument(
        "--composition-radius",
        type=int,
        default=250,
        help=(
            "Radius used to match local GC and target-base density in the "
            "callable background."
        ),
    )
    parser.add_argument(
        "--cluster-bp",
        type=int,
        default=200,
        help="Primary distance for collapsing spatially dependent cases.",
    )
    parser.add_argument(
        "--min-local-controls",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--min-global-controls",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--matches-per-case",
        type=int,
        default=3,
        help=(
            "Maximum nearest controls retained in each case-specific matched "
            "set after hard calipers."
        ),
    )
    parser.add_argument(
        "--min-matches-per-case",
        type=int,
        default=2,
        help=(
            "Minimum controls required in each case-specific matched set."
        ),
    )
    parser.add_argument(
        "--pl-caliper",
        type=float,
        default=0.025,
        help="Maximum absolute control-PL difference within a matched set.",
    )
    parser.add_argument(
        "--log-depth-caliper",
        type=float,
        default=0.30,
        help="Maximum absolute log1p(control-depth) difference.",
    )
    parser.add_argument(
        "--log-support-caliper",
        type=float,
        default=0.30,
        help="Maximum absolute log1p(ME-support) difference.",
    )
    parser.add_argument(
        "--gc-caliper",
        type=float,
        default=0.04,
        help="Maximum absolute local GC-fraction difference.",
    )
    parser.add_argument(
        "--target-caliper",
        type=float,
        default=0.04,
        help="Maximum absolute local target-base-fraction difference.",
    )
    parser.add_argument(
        "--continuous-block-bp",
        type=int,
        default=5000,
        help=(
            "Genomic block size for cluster-robust continuous |delta PL| "
            "models."
        ),
    )
    parser.add_argument(
        "--reservoir-per-stratum",
        type=int,
        default=10000,
        help=(
            "Uniform random reservoir retained per final "
            "contig/strand/callability stratum for genome-wide controls."
        ),
    )
    parser.add_argument(
        "--randomizations",
        type=int,
        default=20000,
    )
    parser.add_argument(
        "--bootstraps",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=170713,
    )
    parser.add_argument(
        "--upstream-bp",
        type=int,
        default=200,
        help="Used only for annotation-sensitivity summaries.",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Retain the temporary SQLite database.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate all inputs and ME coverage, then stop.",
    )
    return parser.parse_args(argv)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def require_directory(path: Path, label: str) -> None:
    if not path.is_dir():
        raise NotADirectoryError(f"{label} is not a directory: {path}")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def safe_float(value: object) -> Optional[float]:
    text = str(value).strip()
    if text.lower() in MISSING:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def stable_seed(seed: int, *parts: object) -> int:
    payload = "|".join([str(seed), *(repr(part) for part in parts)])
    digest = hashlib.blake2b(
        payload.encode("utf-8"),
        digest_size=8,
    ).digest()
    return int.from_bytes(digest, "big") % (2**32 - 1)


def normalize_label(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0 or len(values) == 0:
        return float("nan")
    return float(np.sum(values * weights) / total)


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantile: float,
) -> float:
    if len(values) == 0:
        return float("nan")
    mask = (
        np.isfinite(values)
        & np.isfinite(weights)
        & (weights > 0)
    )
    values = values[mask]
    weights = weights[mask]
    if len(values) == 0:
        return float("nan")
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cutoff = quantile * cumulative[-1]
    index = int(np.searchsorted(cumulative, cutoff, side="left"))
    index = min(index, len(values) - 1)
    return float(values[index])


def effective_sample_size(weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    squared = float(np.sum(weights**2))
    if total <= 0 or squared <= 0:
        return 0.0
    return total * total / squared


def holm_adjust(p_values: Sequence[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    result = np.full(len(p), np.nan)
    finite = np.where(np.isfinite(p))[0]
    if len(finite) == 0:
        return result.tolist()
    order = finite[np.argsort(p[finite])]
    running = 0.0
    m = len(order)
    for rank, index in enumerate(order):
        value = min(1.0, (m - rank) * p[index])
        running = max(running, value)
        result[index] = running
    return result.tolist()


def bh_adjust(p_values: Sequence[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    result = np.full(len(p), np.nan)
    finite = np.where(np.isfinite(p))[0]
    if len(finite) == 0:
        return result.tolist()
    order = finite[np.argsort(p[finite])]
    m = len(order)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = m - reverse_rank + 1
        value = min(1.0, p[index] * m / rank)
        running = min(running, value)
        result[index] = running
    return result.tolist()


def interval_union_length(
    intervals: Iterable[tuple[int, int]],
) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    start, end = ordered[0]
    total = 0
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def configure_logging(output_dir: Path) -> None:
    (output_dir / "logs").mkdir(parents=True, exist_ok=True)
    LOG.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s\t%(levelname)s\t%(message)s"
    )
    file_handler = logging.FileHandler(
        output_dir / "logs" / "analysis.log"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    LOG.handlers.clear()
    LOG.addHandler(file_handler)
    LOG.addHandler(stream_handler)


def collect_upstream_provenance(
    me_inputs: Sequence[tuple[str, str, Path]],
) -> pd.DataFrame:
    """
    Collect, but never invent, chemistry/model provenance. A user-supplied
    label such as '4mC' is not treated as proof that the underlying modBAM and
    entropy alphabet were chemistry-specific.
    """
    rows = []
    for mod_type, sample, directory in me_inputs:
        json_path = directory / "provenance.json"
        log_path = directory / "modkit_entropy.log"
        data: dict[str, object] = {}
        if json_path.is_file():
            try:
                loaded = json.loads(json_path.read_text())
                if isinstance(loaded, dict):
                    data = loaded
            except Exception as error:
                data = {"provenance_json_error": str(error)}

        log_lines = []
        if log_path.is_file():
            try:
                for line in log_path.read_text(
                    errors="replace"
                ).splitlines():
                    lower = line.lower()
                    if any(
                        token in lower
                        for token in (
                            "modkit",
                            "command",
                            "version",
                            "entropy",
                            "filter",
                            "threshold",
                            "motif",
                            "combine-strands",
                            "num-positions",
                            "window-size",
                            "min-coverage",
                        )
                    ):
                        log_lines.append(line.strip())
                    if len(log_lines) >= 20:
                        break
            except OSError as error:
                log_lines = [f"log_read_error={error}"]

        rows.append(
            {
                "mod_type_label": mod_type,
                "control_label": sample,
                "me_directory": str(directory),
                "provenance_json_present": json_path.is_file(),
                "modkit_log_present": log_path.is_file(),
                "dorado_model": data.get("dorado_model", ""),
                "dorado_version": data.get("dorado_version", ""),
                "modkit_version": data.get("modkit_version", ""),
                "modification_code": data.get("modification_code", ""),
                "canonical_base": data.get("canonical_base", ""),
                "entropy_num_positions": data.get(
                    "entropy_num_positions",
                    data.get("num_positions", ""),
                ),
                "entropy_window_size": data.get(
                    "entropy_window_size",
                    data.get("window_size", ""),
                ),
                "entropy_min_coverage": data.get(
                    "entropy_min_coverage",
                    data.get("min_coverage", ""),
                ),
                "combine_strands": data.get("combine_strands", ""),
                "motif": data.get("motif", ""),
                "filter_threshold": data.get("filter_threshold", ""),
                "entropy_command": data.get("entropy_command", ""),
                "selected_log_lines": " || ".join(log_lines),
                "chemistry_specificity_verified_by_script": False,
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# REFERENCE AND OPTIONAL ANNOTATION
# =============================================================================

def read_fai(path: Path) -> OrderedDict[str, int]:
    require_file(path, "FAI")
    result: OrderedDict[str, int] = OrderedDict()
    with path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_number}: malformed FAI row."
                )
            result[fields[0]] = int(fields[1])
    if not result:
        raise ValueError(f"No contigs read from {path}.")
    return result


def read_fasta(path: Path) -> OrderedDict[str, str]:
    require_file(path, "FASTA")
    result: OrderedDict[str, list[str]] = OrderedDict()
    current: Optional[str] = None
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in result:
                    raise ValueError(
                        f"{path}:{line_number}: duplicate contig {current!r}."
                    )
                result[current] = []
            else:
                if current is None:
                    raise ValueError(
                        f"{path}:{line_number}: sequence before FASTA header."
                    )
                result[current].append(line.upper())
    return OrderedDict(
        (chrom, "".join(chunks))
        for chrom, chunks in result.items()
    )


class SequenceReference:
    def __init__(
        self,
        sequences: Mapping[str, str],
        lengths: Mapping[str, int],
    ) -> None:
        self.sequences = dict(sequences)
        self.lengths = dict(lengths)
        self._composition_cache: dict[
            tuple[str, int, str, str, int], tuple[float, float]
        ] = {}
        for chrom, length in lengths.items():
            if chrom not in sequences:
                raise ValueError(f"FASTA missing FAI contig {chrom!r}.")
            if len(sequences[chrom]) != length:
                raise ValueError(
                    f"FASTA/FAI mismatch for {chrom}: "
                    f"{len(sequences[chrom])} != {length}."
                )

    def reference_base(self, key: SiteKey) -> str:
        if key.end - key.start != 1:
            return ""
        base = self.sequences[key.chrom][key.start]
        if key.strand == "-":
            return base.translate(COMPLEMENT)
        return base

    def gc_fraction(self, chrom: str, start: int, end: int) -> float:
        start = max(0, start)
        end = min(self.lengths[chrom], end)
        if end <= start:
            return float("nan")
        seq = self.sequences[chrom][start:end]
        return (seq.count("G") + seq.count("C")) / len(seq)

    def local_composition(
        self,
        key: SiteKey,
        strand_oriented_base: str,
        radius: int,
    ) -> tuple[float, float]:
        cache_key = (
            key.chrom,
            key.start,
            key.strand,
            strand_oriented_base,
            radius,
        )
        cached = self._composition_cache.get(cache_key)
        if cached is not None:
            return cached

        start = max(0, key.start - radius)
        end = min(self.lengths[key.chrom], key.start + radius + 1)
        seq = self.sequences[key.chrom][start:end]
        if not seq:
            result = (float("nan"), float("nan"))
            self._composition_cache[cache_key] = result
            return result

        gc = (seq.count("G") + seq.count("C")) / len(seq)
        target = strand_oriented_base
        if key.strand == "-":
            target = strand_oriented_base.translate(COMPLEMENT)
        target_fraction = (
            seq.count(target) / len(seq)
            if target in DNA_BASES
            else float("nan")
        )
        result = (gc, target_fraction)
        self._composition_cache[cache_key] = result
        return result


def parse_gff_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in text.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            result[key.strip()] = value.strip()
    return result


class AnnotationIndex:
    """
    Broad annotation is used only for diagnostics and a sensitivity analysis.
    No hard-coded organism, prophage or contig coordinates are used.
    """

    def __init__(
        self,
        gff_path: Path,
        operon_path: Path,
        upstream_bp: int,
    ) -> None:
        require_file(gff_path, "GFF")
        require_file(operon_path, "Operon table")
        self.features: dict[
            str, list[tuple[int, int, str, str, str]]
        ] = defaultdict(list)
        gene_features: list[tuple[str, int, int, str, str]] = []

        with gff_path.open() as handle:
            for line_number, raw in enumerate(handle, start=1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                chrom, _, feature_type, start, end, _, strand, _, attrs = fields
                feature_lower = feature_type.lower()
                start0 = int(start) - 1
                end0 = int(end)
                attributes = parse_gff_attributes(attrs)
                locus = (
                    attributes.get("locus_tag")
                    or attributes.get("Name")
                    or attributes.get("ID")
                    or f"{chrom}:{start0}-{end0}"
                )
                if feature_lower in {
                    "rrna", "trna", "ncrna", "tmrna",
                    "srna", "antisense_rna", "riboswitch",
                }:
                    context = "ncRNA"
                elif feature_lower in {"gene", "cds", "pseudogene"}:
                    context = "coding"
                else:
                    continue
                self.features[chrom].append(
                    (start0, end0, context, strand, locus)
                )
                if context == "coding":
                    gene_features.append(
                        (chrom, start0, end0, strand, locus)
                    )

        for chrom in self.features:
            self.features[chrom].sort()

        gene_to_operon = self._read_operons(operon_path)
        groups: dict[str, list[tuple[str, int, int, str, str]]] = defaultdict(list)
        for feature in gene_features:
            locus = feature[4]
            group = gene_to_operon.get(locus, f"gene:{locus}")
            groups[group].append(feature)

        self.upstream: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for members in groups.values():
            chroms = {item[0] for item in members}
            strands = {item[3] for item in members}
            if len(chroms) != 1 or len(strands) != 1:
                continue
            chrom = members[0][0]
            strand = members[0][3]
            if strand == "+":
                tx_start = min(item[1] for item in members)
                self.upstream[chrom].append(
                    (max(0, tx_start - upstream_bp), tx_start)
                )
            elif strand == "-":
                tx_start = max(item[2] for item in members)
                self.upstream[chrom].append(
                    (tx_start, tx_start + upstream_bp)
                )
        for chrom in self.upstream:
            self.upstream[chrom].sort()

    @staticmethod
    def _read_operons(path: Path) -> dict[str, str]:
        frame = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        normalized = {
            normalize_label(column): column
            for column in frame.columns
        }
        op_col = next(
            (
                normalized[name]
                for name in ("operon", "operonid", "idoperon")
                if name in normalized
            ),
            None,
        )
        gene_col = next(
            (
                normalized[name]
                for name in ("idgene", "locustag", "gene", "geneid")
                if name in normalized
            ),
            None,
        )
        if op_col is None or gene_col is None:
            LOG.warning(
                "Operon columns were not recognized; upstream annotation "
                "will use individual genes."
            )
            return {}
        result: dict[str, str] = {}
        for _, row in frame.iterrows():
            locus = str(row[gene_col]).strip()
            operon = str(row[op_col]).strip()
            if locus and operon:
                result[locus] = f"operon:{operon}"
        return result

    def classify(self, key: SiteKey) -> str:
        position = key.start
        overlaps = []
        for start, end, context, _, _ in self.features.get(key.chrom, []):
            if start > position:
                break
            if start <= position < end:
                overlaps.append(context)
        if "ncRNA" in overlaps:
            return "ncRNA"
        if "coding" in overlaps:
            return "coding"
        for start, end in self.upstream.get(key.chrom, []):
            if start > position:
                break
            if start <= position < end:
                return "upstream"
        return "intergenic"


# =============================================================================
# METADATA CASES
# =============================================================================

def parse_mod_cell(
    text: str,
    *,
    path: Path,
    row_number: int,
    column: str,
) -> list[re.Match[str]]:
    text = text.strip()
    if text.lower() in MISSING:
        return []
    tokens = [
        token.strip()
        for token in TOKEN_SPLIT_RE.split(text)
        if token.strip()
    ]
    matches: list[re.Match[str]] = []
    for token in tokens:
        match = SITE_RECORD_RE.fullmatch(token)
        if match is None:
            raise ValueError(
                f"{path}:{row_number}: malformed record in {column!r}: "
                f"{token!r}. Every token must match the complete schema."
            )
        matches.append(match)
    return matches


def read_metadata_cases(
    path: Path,
    mod_labels: set[str],
    contig_lengths: Mapping[str, int],
) -> tuple[
    dict[str, dict[SiteKey, CaseSite]],
    list[str],
    pd.DataFrame,
]:
    require_file(path, "Metadata")
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    columns: list[tuple[str, str, str]] = []
    treatment_order: list[str] = []

    for column in frame.columns:
        column_text = str(column)
        matched = None
        for mod_type in sorted(mod_labels, key=len, reverse=True):
            suffix = "_" + mod_type
            if column_text.endswith(suffix):
                treatment = column_text[:-len(suffix)]
                if treatment:
                    matched = (column_text, treatment, mod_type)
                    break
        if matched is None:
            continue
        columns.append(matched)
        treatment = matched[1]
        if treatment not in treatment_order:
            treatment_order.append(treatment)

    if not columns:
        raise ValueError(
            "No metadata columns matched supplied --me modification labels."
        )

    cases = {
        mod_type: {}
        for mod_type in mod_labels
    }
    qc_rows: list[dict[str, object]] = []

    for row_index, row in frame.iterrows():
        region = str(row.get("region", "")).strip()
        for column, treatment, mod_type in columns:
            for match in parse_mod_cell(
                str(row[column]),
                path=path,
                row_number=row_index + 2,
                column=column,
            ):
                fields = match.groupdict()
                key = SiteKey(
                    chrom=fields["chrom"],
                    start=int(fields["start"]),
                    end=int(fields["end"]),
                    strand=fields["strand"],
                )
                if key.chrom not in contig_lengths:
                    raise ValueError(
                        f"{path}:{row_index + 2}: unknown contig "
                        f"{key.chrom!r}."
                    )
                if not (
                    0 <= key.start < key.end <= contig_lengths[key.chrom]
                ):
                    raise ValueError(
                        f"{path}:{row_index + 2}: out-of-bounds site {key}."
                    )

                reported = float(fields["reported_delta"])
                control_pl = (
                    float(fields["control_pl"])
                    if fields["control_pl"] is not None
                    else None
                )
                treatment_pl = (
                    float(fields["treatment_pl"])
                    if fields["treatment_pl"] is not None
                    else None
                )
                for label, value in (
                    ("control", control_pl),
                    ("treatment", treatment_pl),
                ):
                    if value is not None and not 0 <= value <= 1:
                        raise ValueError(
                            f"{path}:{row_index + 2}: {label} PL {value} "
                            "is outside 0..1."
                        )

                case = cases[mod_type].setdefault(
                    key,
                    CaseSite(mod_type=mod_type, key=key),
                )
                case.treatments.add(treatment)
                case.reported_deltas[treatment].append(reported)
                if control_pl is not None and treatment_pl is not None:
                    case.control_pl[treatment].append(control_pl)
                    case.treatment_pl[treatment].append(treatment_pl)
                    calculated = treatment_pl - control_pl
                    sign_relation = (
                        "same"
                        if math.copysign(1.0, calculated)
                        == math.copysign(1.0, reported)
                        else "opposite"
                    )
                    discrepancy_same = abs(reported - calculated)
                    discrepancy_opposite = abs(reported + calculated)
                else:
                    calculated = float("nan")
                    sign_relation = "unavailable"
                    discrepancy_same = float("nan")
                    discrepancy_opposite = float("nan")
                case.source_regions.add(region)

                qc_rows.append(
                    {
                        "row": row_index + 2,
                        "column": column,
                        "mod_type": mod_type,
                        "treatment": treatment,
                        "chrom": key.chrom,
                        "start": key.start,
                        "end": key.end,
                        "strand": key.strand,
                        "reported_delta": reported,
                        "control_pl": control_pl,
                        "treatment_pl": treatment_pl,
                        "calculated_treatment_minus_control": calculated,
                        "reported_sign_relation": sign_relation,
                        "abs_error_same_sign": discrepancy_same,
                        "abs_error_opposite_sign": discrepancy_opposite,
                    }
                )

    return cases, treatment_order, pd.DataFrame(qc_rows)


# =============================================================================
# DMR DISCOVERY AND SQLITE RISK SET
# =============================================================================

def best_treatment_match(
    path_label: str,
    treatments: Sequence[str],
) -> Optional[str]:
    normalized = normalize_label(path_label)
    direct = {
        normalize_label(treatment): treatment
        for treatment in treatments
    }
    if normalized in direct:
        return direct[normalized]

    # Common visual O/0 confusion is handled generically.
    variants = {
        normalized,
        normalized.replace("0", "o"),
        normalized.replace("o", "0"),
    }
    for variant in variants:
        for key, treatment in direct.items():
            if variant == key or variant == key.replace("0", "o"):
                return treatment

    scores = [
        (
            difflib.SequenceMatcher(
                None,
                normalized,
                normalize_label(treatment),
            ).ratio(),
            treatment,
        )
        for treatment in treatments
    ]
    if not scores:
        return None
    scores.sort(reverse=True)
    if scores[0][0] >= 0.72:
        if len(scores) == 1 or scores[0][0] - scores[1][0] >= 0.10:
            return scores[0][1]
    return None


def discover_dmr_files(
    root: Path,
    mod_labels: set[str],
    treatments: Sequence[str],
) -> dict[tuple[str, str], Path]:
    require_directory(root, "DMR root")
    result: dict[tuple[str, str], Path] = {}

    for path in sorted(root.rglob("control_treatment.bed")):
        parts = list(path.parts[:-1])
        mod_type = next(
            (
                label
                for label in mod_labels
                if normalize_label(label)
                in {normalize_label(part) for part in parts}
            ),
            None,
        )
        if mod_type is None:
            continue

        treatment = None
        for part in reversed(parts):
            treatment = best_treatment_match(part, treatments)
            if treatment is not None:
                break
        if treatment is None:
            LOG.warning(
                "Could not infer treatment for DMR file %s; ignored.",
                path,
            )
            continue

        key = (mod_type, treatment)
        if key in result:
            raise ValueError(
                f"Multiple DMR files map to {key}: "
                f"{result[key]} and {path}."
            )
        result[key] = path

    missing = [
        f"{mod}/{treatment}"
        for mod in sorted(mod_labels)
        for treatment in treatments
        if (mod, treatment) not in result
    ]
    if missing:
        LOG.warning(
            "No DMR file was found for: %s",
            ", ".join(missing),
        )
    return result


def sqlite_table_name(mod_type: str) -> str:
    digest = hashlib.blake2b(
        mod_type.encode("utf-8"),
        digest_size=5,
    ).hexdigest()
    return f"sites_{digest}"


def create_site_table(
    connection: sqlite3.Connection,
    mod_type: str,
    treatment_count: int,
) -> str:
    table = sqlite_table_name(mod_type)
    dynamic = []
    for index in range(treatment_count):
        dynamic.extend(
            [
                f"pl_{index} REAL",
                f"depth_{index} REAL",
                f"post_pl_{index} REAL",
                f"post_depth_{index} REAL",
            ]
        )
    connection.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {table} (
            chrom TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT NOT NULL,
            name TEXT NOT NULL,
            mask INTEGER NOT NULL DEFAULT 0,
            {", ".join(dynamic)},
            PRIMARY KEY (chrom, start, end, strand)
        )
        """
    )
    connection.execute(
        f"CREATE INDEX IF NOT EXISTS idx_{table}_range "
        f"ON {table}(chrom, strand, start)"
    )
    connection.execute(
        f"CREATE INDEX IF NOT EXISTS idx_{table}_mask "
        f"ON {table}(mask)"
    )
    return table


def dmr_header(path: Path) -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row and any(field.strip() for field in row):
                return [field.strip() for field in row]
    raise ValueError(f"{path}: empty DMR file.")


def stream_dmr_into_sqlite(
    connection: sqlite3.Connection,
    dmr: DMRFile,
    treatment_count: int,
    chunk_size: int = 25000,
) -> dict[str, object]:
    required = {
        "chrom", "start", "end", "name", "strand",
        "control_total", "treatment_total",
        "control_pct_modified", "treatment_pct_modified",
    }
    header = dmr_header(dmr.path)
    missing = required - set(header)
    if missing:
        raise ValueError(
            f"{dmr.path}: missing DMR columns {sorted(missing)}."
        )

    table = sqlite_table_name(dmr.mod_type)
    columns = [
        "chrom", "start", "end", "strand", "name", "mask",
        f"pl_{dmr.treatment_index}",
        f"depth_{dmr.treatment_index}",
        f"post_pl_{dmr.treatment_index}",
        f"post_depth_{dmr.treatment_index}",
    ]
    placeholders = ",".join("?" for _ in columns)
    updates = [
        "mask = mask | excluded.mask",
        "name = CASE WHEN name = excluded.name THEN name ELSE name END",
        f"pl_{dmr.treatment_index} = excluded.pl_{dmr.treatment_index}",
        f"depth_{dmr.treatment_index} = excluded.depth_{dmr.treatment_index}",
        f"post_pl_{dmr.treatment_index} = excluded.post_pl_{dmr.treatment_index}",
        f"post_depth_{dmr.treatment_index} = excluded.post_depth_{dmr.treatment_index}",
    ]
    sql = (
        f"INSERT INTO {table} ({','.join(columns)}) "
        f"VALUES ({placeholders}) "
        f"ON CONFLICT(chrom,start,end,strand) DO UPDATE SET "
        f"{', '.join(updates)}"
    )

    rows = 0
    valid_pl = 0
    name_counts: Counter[str] = Counter()
    batch: list[tuple[object, ...]] = []

    with open_text(dmr.path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=2):
            if not row or all(
                str(value).strip() == ""
                for value in row.values()
            ):
                continue
            chrom = str(row["chrom"]).strip()
            start = int(row["start"])
            end = int(row["end"])
            strand = str(row["strand"]).strip() or "."
            name = str(row["name"]).strip().upper()
            control_pl = safe_float(row["control_pct_modified"])
            control_depth = safe_float(row["control_total"])
            post_pl = safe_float(row["treatment_pct_modified"])
            post_depth = safe_float(row["treatment_total"])

            if control_pl is not None and not 0 <= control_pl <= 1:
                raise ValueError(
                    f"{dmr.path}:{row_number}: control PL outside 0..1."
                )
            if post_pl is not None and not 0 <= post_pl <= 1:
                raise ValueError(
                    f"{dmr.path}:{row_number}: treatment PL outside 0..1."
                )
            if control_pl is not None:
                valid_pl += 1
            if name in DNA_BASES:
                name_counts[name] += 1

            batch.append(
                (
                    chrom,
                    start,
                    end,
                    strand,
                    name,
                    dmr.bit,
                    control_pl,
                    control_depth,
                    post_pl,
                    post_depth,
                )
            )
            rows += 1
            if len(batch) >= chunk_size:
                connection.executemany(sql, batch)
                connection.commit()
                batch.clear()

    if batch:
        connection.executemany(sql, batch)
        connection.commit()

    dominant_base = (
        name_counts.most_common(1)[0][0]
        if name_counts
        else ""
    )
    return {
        "mod_type": dmr.mod_type,
        "treatment": dmr.treatment,
        "path": str(dmr.path),
        "rows": rows,
        "rows_with_control_pl": valid_pl,
        "dominant_name_base": dominant_base,
    }


def fetch_site_record(
    connection: sqlite3.Connection,
    mod_type: str,
    key: SiteKey,
    treatment_count: int,
) -> Optional[SiteRecord]:
    table = sqlite_table_name(mod_type)
    dynamic = []
    for index in range(treatment_count):
        dynamic.extend(
            [
                f"pl_{index}",
                f"depth_{index}",
                f"post_pl_{index}",
                f"post_depth_{index}",
            ]
        )
    row = connection.execute(
        f"""
        SELECT chrom,start,end,strand,name,mask,{",".join(dynamic)}
        FROM {table}
        WHERE chrom=? AND start=? AND end=? AND strand=?
        """,
        (key.chrom, key.start, key.end, key.strand),
    ).fetchone()
    if row is None:
        return None
    return row_to_site_record(row, treatment_count)


def row_to_site_record(
    row: Sequence[object],
    treatment_count: int,
    sampling_weight: float = 1.0,
) -> SiteRecord:
    chrom, start, end, strand, name, mask = row[:6]
    dynamic = row[6:]
    pls = []
    depths = []
    post_pls = []
    post_depths = []
    for index in range(treatment_count):
        offset = index * 4
        pls.append(safe_float(dynamic[offset]))
        depths.append(safe_float(dynamic[offset + 1]))
        post_pls.append(safe_float(dynamic[offset + 2]))
        post_depths.append(safe_float(dynamic[offset + 3]))
    return SiteRecord(
        key=SiteKey(str(chrom), int(start), int(end), str(strand)),
        name=str(name),
        mask=int(mask),
        pls=tuple(pls),
        depths=tuple(depths),
        post_pls=tuple(post_pls),
        post_depths=tuple(post_depths),
        sampling_weight=sampling_weight,
    )


def relevant_average(
    values: Sequence[Optional[float]],
    required_indices: Sequence[int],
) -> Optional[float]:
    selected = [
        values[index]
        for index in required_indices
        if index < len(values)
        and values[index] is not None
        and math.isfinite(float(values[index]))
    ]
    if not selected:
        return None
    return float(np.mean(selected))


# =============================================================================
# MODKIT ME INPUTS
# =============================================================================

def parse_windows_bedgraph(
    mod_type: str,
    sample: str,
    path: Path,
    contig_lengths: Mapping[str, int],
) -> tuple[list[WindowTrack], list[dict[str, object]]]:
    """
    Parse genuine Modkit entropy windows.

    Accepted layouts are:

        chrom start end entropy support
        chrom start end entropy strand support
        chrom start end strand entropy support
        chrom start end entropy support strand

    The second layout is the current study's Modkit output. Five-column input
    is explicitly unstranded. Duplicate intervals are rejected because entropy
    cannot be pooled by averaging scalar entropy estimates.
    """
    require_file(path, "windows.bedgraph")
    grouped: dict[
        tuple[str, str],
        list[tuple[int, int, float, float]]
    ] = defaultdict(list)
    seen: set[tuple[str, int, int, str]] = set()
    total_rows = 0
    finite_rows = 0
    zero_support_rows = 0
    malformed = 0
    detected_layouts: Counter[str] = Counter()

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()

            if len(fields) >= 3:
                f0, f1, f2 = (
                    fields[0].strip().lower(),
                    fields[1].strip().lower(),
                    fields[2].strip().lower(),
                )
                if (
                    f0 in {"chrom", "chromosome", "contig"}
                    and f1 in {"start", "chromstart"}
                    and f2 in {"end", "chromend", "stop"}
                ):
                    continue

            if len(fields) < 5:
                malformed += 1
                raise ValueError(
                    f"{path}:{line_number}: expected at least five columns."
                )

            chrom = fields[0].strip()
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: start/end are not integers."
                ) from error

            strand = "."
            entropy: Optional[float]
            support: Optional[float]

            if (
                len(fields) >= 6
                and fields[4].strip() in {"+", "-", "."}
            ):
                # Modkit layout used in this study:
                # chrom start end entropy strand read_support
                entropy = safe_float(fields[3])
                strand = fields[4].strip()
                support = safe_float(fields[5])
                detected_layouts["entropy_strand_support"] += 1
            elif (
                len(fields) >= 6
                and fields[3].strip() in {"+", "-", "."}
            ):
                # chrom start end strand entropy read_support
                strand = fields[3].strip()
                entropy = safe_float(fields[4])
                support = safe_float(fields[5])
                detected_layouts["strand_entropy_support"] += 1
            elif (
                len(fields) >= 6
                and fields[5].strip() in {"+", "-", "."}
            ):
                # chrom start end entropy read_support strand
                entropy = safe_float(fields[3])
                support = safe_float(fields[4])
                strand = fields[5].strip()
                detected_layouts["entropy_support_strand"] += 1
            else:
                # Standard unstranded five-column bedGraph.
                entropy = safe_float(fields[3])
                support = safe_float(fields[4])
                detected_layouts["entropy_support_unstranded"] += 1

            total_rows += 1
            if chrom not in contig_lengths:
                raise ValueError(
                    f"{path}:{line_number}: unknown contig {chrom!r}."
                )
            if not (0 <= start < end <= contig_lengths[chrom]):
                raise ValueError(
                    f"{path}:{line_number}: invalid interval "
                    f"{chrom}:{start}-{end}."
                )
            if entropy is None or support is None:
                continue
            if not 0 <= entropy <= 1:
                raise ValueError(
                    f"{path}:{line_number}: entropy outside 0..1."
                )
            if support < 0:
                raise ValueError(
                    f"{path}:{line_number}: negative support."
                )
            if support == 0:
                zero_support_rows += 1

            duplicate_key = (chrom, start, end, strand)
            if duplicate_key in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate entropy interval "
                    f"{duplicate_key}. Duplicates are rejected because entropy "
                    "cannot be combined by averaging."
                )
            seen.add(duplicate_key)
            grouped[(chrom, strand)].append(
                (start, end, entropy, support)
            )
            finite_rows += 1

    tracks: list[WindowTrack] = []
    qc_rows: list[dict[str, object]] = []
    layout_summary = ",".join(
        f"{name}:{count}"
        for name, count in sorted(detected_layouts.items())
    )

    for (chrom, strand), rows in grouped.items():
        rows.sort(key=lambda item: (item[0], item[1]))
        starts = np.asarray([row[0] for row in rows], dtype=np.int64)
        ends = np.asarray([row[1] for row in rows], dtype=np.int64)
        entropy_values = np.asarray(
            [row[2] for row in rows],
            dtype=float,
        )
        support_values = np.asarray(
            [row[3] for row in rows],
            dtype=float,
        )
        midpoints = (starts + ends) / 2.0
        widths = ends - starts
        steps = np.diff(midpoints)

        tracks.append(
            WindowTrack(
                sample=sample,
                mod_type=mod_type,
                path=path,
                stranded=(strand != "."),
                chrom=chrom,
                strand=strand,
                starts=starts,
                ends=ends,
                midpoints=midpoints,
                entropy=entropy_values,
                support=support_values,
                median_width=float(np.median(widths)),
                median_step=(
                    float(np.median(steps[steps > 0]))
                    if np.any(steps > 0)
                    else float(np.median(widths))
                ),
                max_width=int(np.max(widths)),
            )
        )

        finite_coverage = interval_union_length(
            zip(starts.tolist(), ends.tolist())
        )
        qc_rows.append(
            {
                "mod_type": mod_type,
                "sample": sample,
                "file": str(path),
                "chrom": chrom,
                "strand": strand,
                "format": "windows.bedgraph",
                "detected_column_layouts": layout_summary,
                "total_file_rows": total_rows,
                "finite_file_rows": finite_rows,
                "finite_rows": len(rows),
                "malformed_rows": malformed,
                "zero_support_rows_total": zero_support_rows,
                "finite_coverage_bp": finite_coverage,
                "contig_length_bp": contig_lengths[chrom],
                "finite_coverage_fraction": (
                    finite_coverage / contig_lengths[chrom]
                ),
                "median_window_width_bp": float(np.median(widths)),
                "median_window_step_bp": (
                    float(np.median(steps[steps > 0]))
                    if np.any(steps > 0)
                    else float("nan")
                ),
                "median_support": float(np.median(support_values)),
                "min_support": float(np.min(support_values)),
                "max_support": float(np.max(support_values)),
            }
        )

    if not tracks:
        raise ValueError(f"No finite entropy windows read from {path}.")
    return tracks, qc_rows


def parse_regions_qc(
    mod_type: str,
    sample: str,
    path: Path,
    contig_lengths: Mapping[str, int],
) -> list[dict[str, object]]:
    require_file(path, "regions.bed")
    requested: dict[str, list[tuple[int, int]]] = defaultdict(list)
    finite: dict[str, list[tuple[int, int]]] = defaultdict(list)
    rows_by_contig: Counter[str] = Counter()
    finite_by_contig: Counter[str] = Counter()
    failed_by_contig: Counter[str] = Counter()
    successful_by_contig: Counter[str] = Counter()
    seen: set[tuple[str, int, int, str, str]] = set()

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")

            if len(fields) >= 3:
                f0, f1, f2 = (
                    fields[0].strip().lower(),
                    fields[1].strip().lower(),
                    fields[2].strip().lower(),
                )
                if (
                    f0 in {"chrom", "chromosome", "contig"}
                    and f1 in {"start", "chromstart"}
                    and f2 in {"end", "chromend", "stop"}
                ):
                    continue

            if len(fields) < 14:
                raise ValueError(
                    f"{path}:{line_number}: regions.bed must have 14 columns."
                )
            chrom = fields[0].strip()
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3].strip()
            mean_entropy = safe_float(fields[4])
            strand = fields[5].strip() or "."
            successful = int(float(fields[12]))
            failed = int(float(fields[13]))

            if chrom not in contig_lengths:
                raise ValueError(
                    f"{path}:{line_number}: unknown contig {chrom!r}."
                )
            if not (0 <= start < end <= contig_lengths[chrom]):
                raise ValueError(
                    f"{path}:{line_number}: invalid region interval."
                )
            key = (chrom, start, end, strand, name)
            if key in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate regions.bed row {key}."
                )
            seen.add(key)

            requested[chrom].append((start, end))
            rows_by_contig[chrom] += 1
            successful_by_contig[chrom] += successful
            failed_by_contig[chrom] += failed
            if mean_entropy is not None:
                finite[chrom].append((start, end))
                finite_by_contig[chrom] += 1

    result = []
    for chrom, length in contig_lengths.items():
        requested_bp = interval_union_length(requested.get(chrom, []))
        finite_bp = interval_union_length(finite.get(chrom, []))
        result.append(
            {
                "mod_type": mod_type,
                "sample": sample,
                "file": str(path),
                "chrom": chrom,
                "format": "regions.bed",
                "requested_rows": rows_by_contig[chrom],
                "finite_entropy_rows": finite_by_contig[chrom],
                "requested_coverage_bp": requested_bp,
                "finite_entropy_coverage_bp": finite_bp,
                "contig_length_bp": length,
                "requested_coverage_fraction": requested_bp / length,
                "finite_entropy_coverage_fraction": finite_bp / length,
                "successful_window_count": successful_by_contig[chrom],
                "failed_window_count": failed_by_contig[chrom],
                "successful_fraction": (
                    successful_by_contig[chrom]
                    / (
                        successful_by_contig[chrom]
                        + failed_by_contig[chrom]
                    )
                    if successful_by_contig[chrom]
                    + failed_by_contig[chrom]
                    > 0
                    else float("nan")
                ),
            }
        )
    return result


class MEIndex:
    def __init__(
        self,
        tracks: Sequence[WindowTrack],
    ) -> None:
        self.tracks: dict[
            tuple[str, str, str, str], WindowTrack
        ] = {}
        self.samples_by_mod: dict[str, list[str]] = defaultdict(list)
        self.strand_modes: dict[str, set[str]] = defaultdict(set)
        for track in tracks:
            key = (
                track.mod_type,
                track.sample,
                track.chrom,
                track.strand,
            )
            if key in self.tracks:
                raise ValueError(f"Duplicate ME track key: {key}")
            self.tracks[key] = track
            if track.sample not in self.samples_by_mod[track.mod_type]:
                self.samples_by_mod[track.mod_type].append(track.sample)
            self.strand_modes[track.mod_type].add(track.strand)

    def track_for(
        self,
        mod_type: str,
        sample: str,
        key: SiteKey,
    ) -> tuple[Optional[WindowTrack], str]:
        exact = self.tracks.get(
            (mod_type, sample, key.chrom, key.strand)
        )
        if exact is not None:
            return exact, "same_strand"
        unstranded = self.tracks.get(
            (mod_type, sample, key.chrom, ".")
        )
        if unstranded is not None:
            return unstranded, "unstranded"
        return None, "missing"

    @staticmethod
    def _overlap_indices(
        track: WindowTrack,
        position: int,
    ) -> np.ndarray:
        left = int(np.searchsorted(
            track.starts,
            position - track.max_width,
            side="left",
        ))
        right = int(np.searchsorted(
            track.starts,
            position,
            side="right",
        ))
        if right <= left:
            return np.array([], dtype=int)
        candidate = np.arange(left, right, dtype=int)
        return candidate[track.ends[left:right] > position]

    @staticmethod
    def _midpoint_indices(
        track: WindowTrack,
        position: int,
        radius: int,
    ) -> np.ndarray:
        left = int(np.searchsorted(
            track.midpoints,
            position - radius,
            side="left",
        ))
        right = int(np.searchsorted(
            track.midpoints,
            position + radius,
            side="right",
        ))
        return np.arange(left, right, dtype=int)

    def score(
        self,
        mod_type: str,
        key: SiteKey,
        min_support: int,
    ) -> Optional[SiteScore]:
        sample_scores = []
        strand_modes = []

        for sample in self.samples_by_mod.get(mod_type, []):
            track, strand_mode = self.track_for(
                mod_type,
                sample,
                key,
            )
            if track is None:
                continue
            overlap = self._overlap_indices(track, key.position)
            overlap = overlap[
                track.support[overlap] >= min_support
            ]
            if len(overlap) == 0:
                continue

            site_values = track.entropy[overlap]
            site_weights = track.support[overlap]
            site_me = weighted_quantile(
                site_values,
                site_weights,
                0.50,
            )
            site_q75 = weighted_quantile(
                site_values,
                site_weights,
                0.75,
            )
            site_support = weighted_mean(
                track.support[overlap],
                np.ones(len(overlap), dtype=float),
            )

            domain_values = {}
            for radius in (250, 500, 1000):
                indices = self._midpoint_indices(
                    track,
                    key.position,
                    radius,
                )
                indices = indices[
                    track.support[indices] >= min_support
                ]
                domain_values[radius] = (
                    weighted_quantile(
                        track.entropy[indices],
                        track.support[indices],
                        0.50,
                    )
                    if len(indices) > 0
                    else float("nan")
                )

            # Direct focal-peak contrast requested by the study:
            # compare the entropy windows physically overlapping the base with
            # non-overlapping entropy windows whose midpoints fall within
            # +/-500 bp of that same base.
            overlap_set = set(overlap.tolist())
            context_500_indices = self._midpoint_indices(
                track,
                key.position,
                500,
            )
            context_500_indices = context_500_indices[
                track.support[context_500_indices] >= min_support
            ]
            context_500_indices = np.asarray(
                [
                    index
                    for index in context_500_indices
                    if index not in overlap_set
                ],
                dtype=int,
            )
            context_500_median = (
                weighted_quantile(
                    track.entropy[context_500_indices],
                    track.support[context_500_indices],
                    0.50,
                )
                if len(context_500_indices) >= 3
                else float("nan")
            )
            local_prominence_500 = (
                site_me - context_500_median
                if math.isfinite(context_500_median)
                else float("nan")
            )
            local_percentile_500 = (
                float(np.sum(
                    track.support[context_500_indices][
                        track.entropy[context_500_indices] <= site_me
                    ]
                ))
                / float(np.sum(track.support[context_500_indices]))
                if (
                    len(context_500_indices) >= 3
                    and float(np.sum(track.support[context_500_indices])) > 0
                )
                else float("nan")
            )

            # Retain the historical +/-1000-bp descriptive contrast for
            # backward compatibility with Versions 3/4.
            local_indices = self._midpoint_indices(
                track,
                key.position,
                1000,
            )
            local_indices = local_indices[
                track.support[local_indices] >= min_support
            ]
            baseline_indices = np.asarray(
                [
                    index
                    for index in local_indices
                    if index not in overlap_set
                ],
                dtype=int,
            )
            local_median = (
                weighted_quantile(
                    track.entropy[baseline_indices],
                    track.support[baseline_indices],
                    0.50,
                )
                if len(baseline_indices) >= 3
                else float("nan")
            )
            prominence = (
                site_me - local_median
                if math.isfinite(local_median)
                else float("nan")
            )

            sample_scores.append(
                {
                    "site_me": site_me,
                    "site_q75": site_q75,
                    "site_support": site_support,
                    "overlap_windows": len(overlap),
                    "domain_250": domain_values[250],
                    "domain_500": domain_values[500],
                    "domain_1000": domain_values[1000],
                    "context_500_excluding_site": context_500_median,
                    "local_prominence_500": local_prominence_500,
                    "local_percentile_500": local_percentile_500,
                    "local_context_windows_500": len(context_500_indices),
                    "local_window_median": local_median,
                    "local_window_prominence": prominence,
                }
            )
            strand_modes.append(strand_mode)

        if not sample_scores:
            return None

        def median_metric(name: str) -> float:
            values = np.asarray(
                [row[name] for row in sample_scores],
                dtype=float,
            )
            finite = values[np.isfinite(values)]
            return (
                float(np.median(finite))
                if len(finite)
                else float("nan")
            )

        strand_mode = (
            "same_strand"
            if all(mode == "same_strand" for mode in strand_modes)
            else "unstranded"
        )
        return SiteScore(
            mod_type=mod_type,
            key=key,
            strand_mode=strand_mode,
            site_me=median_metric("site_me"),
            site_q75=median_metric("site_q75"),
            site_support=median_metric("site_support"),
            overlap_windows=int(round(
                median_metric("overlap_windows")
            )),
            domain_250=median_metric("domain_250"),
            domain_500=median_metric("domain_500"),
            domain_1000=median_metric("domain_1000"),
            context_500_excluding_site=median_metric(
                "context_500_excluding_site"
            ),
            local_prominence_500=median_metric(
                "local_prominence_500"
            ),
            local_percentile_500=median_metric(
                "local_percentile_500"
            ),
            local_context_windows_500=int(round(
                median_metric("local_context_windows_500")
            )),
            local_window_median=median_metric("local_window_median"),
            local_window_prominence=median_metric(
                "local_window_prominence"
            ),
            control_samples=len(sample_scores),
        )


# =============================================================================
# RANDOM RISK-SET RESERVOIRS
# =============================================================================

class Reservoir:
    def __init__(self, maximum: int, seed: int) -> None:
        self.maximum = maximum
        self.rng = random.Random(seed)
        self.seen = 0
        self.items: list[SiteRecord] = []

    def add(self, item: SiteRecord) -> None:
        self.seen += 1
        if len(self.items) < self.maximum:
            self.items.append(item)
            return
        index = self.rng.randrange(self.seen)
        if index < self.maximum:
            self.items[index] = item

    def finalize(self) -> None:
        weight = (
            self.seen / len(self.items)
            if self.items
            else 0.0
        )
        for item in self.items:
            item.sampling_weight = weight


def build_global_reservoirs(
    connection: sqlite3.Connection,
    mod_type: str,
    treatment_count: int,
    case_keys: set[SiteKey],
    maximum: int,
    seed: int,
) -> tuple[
    dict[tuple[str, str, int], Reservoir],
    pd.DataFrame,
]:
    table = sqlite_table_name(mod_type)
    dynamic = []
    for index in range(treatment_count):
        dynamic.extend(
            [
                f"pl_{index}",
                f"depth_{index}",
                f"post_pl_{index}",
                f"post_depth_{index}",
            ]
        )
    query = (
        f"SELECT chrom,start,end,strand,name,mask,{','.join(dynamic)} "
        f"FROM {table}"
    )

    reservoirs: dict[tuple[str, str, int], Reservoir] = {}
    cursor = connection.execute(query)
    for row in cursor:
        record = row_to_site_record(row, treatment_count)
        if record.key in case_keys:
            continue
        stratum = (
            record.key.chrom,
            record.key.strand,
            record.mask,
        )
        reservoir = reservoirs.setdefault(
            stratum,
            Reservoir(
                maximum,
                stable_seed(seed, mod_type, stratum),
            ),
        )
        reservoir.add(record)

    qc_rows = []
    for stratum, reservoir in reservoirs.items():
        reservoir.finalize()
        qc_rows.append(
            {
                "mod_type": mod_type,
                "chrom": stratum[0],
                "strand": stratum[1],
                "callability_mask": stratum[2],
                "eligible_background_rows": reservoir.seen,
                "reservoir_rows": len(reservoir.items),
                "sampling_fraction": (
                    len(reservoir.items) / reservoir.seen
                    if reservoir.seen
                    else 0.0
                ),
                "inverse_sampling_weight": (
                    reservoir.seen / len(reservoir.items)
                    if reservoir.items
                    else float("nan")
                ),
            }
        )
    return reservoirs, pd.DataFrame(qc_rows)


# =============================================================================
# CASE/CONTROL POOLS
# =============================================================================

def required_indices(
    case: CaseSite,
    treatment_to_index: Mapping[str, int],
) -> list[int]:
    return sorted(
        treatment_to_index[treatment]
        for treatment in case.treatments
        if treatment in treatment_to_index
    )


def case_mask(
    case: CaseSite,
    treatment_to_index: Mapping[str, int],
) -> int:
    mask = 0
    for treatment in case.treatments:
        if treatment in treatment_to_index:
            mask |= 1 << treatment_to_index[treatment]
    return mask


def case_covariates(
    case: CaseSite,
    record: SiteRecord,
    treatment_to_index: Mapping[str, int],
) -> tuple[Optional[float], Optional[float]]:
    indices = required_indices(case, treatment_to_index)
    metadata_values = [
        value
        for treatment in case.treatments
        for value in case.control_pl.get(treatment, [])
        if value is not None
    ]
    control_pl = (
        float(np.mean(metadata_values))
        if metadata_values
        else relevant_average(record.pls, indices)
    )
    depth = relevant_average(record.depths, indices)
    return control_pl, depth


def strict_match_weights(
    case_pl: float,
    case_depth: float,
    case_support: float,
    case_gc: float,
    case_target_fraction: float,
    controls_pl: np.ndarray,
    controls_depth: np.ndarray,
    controls_support: np.ndarray,
    controls_gc: np.ndarray,
    controls_target_fraction: np.ndarray,
    base_weights: np.ndarray,
    *,
    matches_per_case: int,
    min_matches_per_case: int,
    pl_caliper: float,
    log_depth_caliper: float,
    log_support_caliper: float,
    gc_caliper: float,
    target_caliper: float,
) -> tuple[np.ndarray, float]:
    """
    Case-specific nearest-neighbour matching with replacement across cases.

    Controls must pass every hard caliper. Within the admissible set, a small
    balanced subset is chosen greedily so that the mean standardized
    covariate difference of the matched set remains close to zero. Outcomes
    are never used for selection.

    The returned array has equal weight on selected controls and zero
    elsewhere. The second return value is the maximum absolute standardized
    mean mismatch within that case-specific set.
    """
    n = len(controls_pl)
    empty = np.zeros(n, dtype=float)
    if n == 0:
        return empty, float("nan")

    case_log_depth = math.log1p(case_depth)
    case_log_support = math.log1p(case_support)
    control_log_depth = np.log1p(controls_depth)
    control_log_support = np.log1p(controls_support)

    raw_differences = np.column_stack(
        [
            controls_pl - case_pl,
            control_log_depth - case_log_depth,
            control_log_support - case_log_support,
            controls_gc - case_gc,
            controls_target_fraction - case_target_fraction,
        ]
    )
    calipers = np.asarray(
        [
            pl_caliper,
            log_depth_caliper,
            log_support_caliper,
            gc_caliper,
            target_caliper,
        ],
        dtype=float,
    )

    finite = np.all(np.isfinite(raw_differences), axis=1)
    admissible = finite & np.all(
        np.abs(raw_differences) <= calipers,
        axis=1,
    )
    candidate_indices = np.flatnonzero(admissible)
    if len(candidate_indices) < min_matches_per_case:
        return empty, float("nan")

    standardized = raw_differences[candidate_indices] / calipers
    individual_distance = np.sum(standardized ** 2, axis=1)

    target_count = min(matches_per_case, len(candidate_indices))
    selected_local: list[int] = []
    selected_mask = np.zeros(len(candidate_indices), dtype=bool)
    running_sum = np.zeros(standardized.shape[1], dtype=float)

    # Greedily minimize the mean covariate mismatch of the matched set while
    # retaining a small penalty for individually distant controls. The tiny
    # sampling-weight term only breaks exact ties reproducibly.
    safe_base = np.asarray(base_weights, dtype=float)[candidate_indices]
    safe_base = np.where(
        np.isfinite(safe_base) & (safe_base > 0),
        safe_base,
        1.0,
    )
    tie_term = -1e-12 * np.log(safe_base)

    for step in range(target_count):
        available = ~selected_mask
        trial_mean = (
            running_sum[None, :]
            + standardized
        ) / float(step + 1)
        objective = (
            np.sum(trial_mean ** 2, axis=1)
            + 0.05 * individual_distance
            + tie_term
        )
        objective[~available] = np.inf
        chosen = int(np.argmin(objective))
        if not math.isfinite(float(objective[chosen])):
            break
        selected_mask[chosen] = True
        selected_local.append(chosen)
        running_sum += standardized[chosen]

    if len(selected_local) < min_matches_per_case:
        return empty, float("nan")

    selected_global = candidate_indices[
        np.asarray(selected_local, dtype=int)
    ]
    weights = np.zeros(n, dtype=float)
    weights[selected_global] = 1.0
    mean_standardized_difference = (
        running_sum / float(len(selected_local))
    )
    maximum_mismatch = float(
        np.max(np.abs(mean_standardized_difference))
    )
    return weights, maximum_mismatch


def query_local_controls(
    connection: sqlite3.Connection,
    mod_type: str,
    case: CaseSite,
    treatment_count: int,
    treatment_to_index: Mapping[str, int],
    radius: int,
    all_case_keys: set[SiteKey],
) -> list[SiteRecord]:
    table = sqlite_table_name(mod_type)
    required = case_mask(case, treatment_to_index)
    dynamic = []
    for index in range(treatment_count):
        dynamic.extend(
            [
                f"pl_{index}",
                f"depth_{index}",
                f"post_pl_{index}",
                f"post_depth_{index}",
            ]
        )

    rows = connection.execute(
        f"""
        SELECT chrom,start,end,strand,name,mask,{",".join(dynamic)}
        FROM {table}
        WHERE chrom=?
          AND strand=?
          AND start BETWEEN ? AND ?
          AND (mask & ?) = ?
        """,
        (
            case.key.chrom,
            case.key.strand,
            max(0, case.key.start - radius),
            case.key.start + radius,
            required,
            required,
        ),
    )
    controls = []
    for row in rows:
        record = row_to_site_record(row, treatment_count)
        if record.key in all_case_keys:
            continue
        controls.append(record)
    return controls


def eligible_global_controls(
    reservoirs: Mapping[tuple[str, str, int], Reservoir],
    case: CaseSite,
    treatment_to_index: Mapping[str, int],
) -> list[SiteRecord]:
    required = case_mask(case, treatment_to_index)
    controls = []
    for (chrom, strand, mask), reservoir in reservoirs.items():
        if chrom != case.key.chrom or strand != case.key.strand:
            continue
        if (mask & required) != required:
            continue
        controls.extend(reservoir.items)
    return controls


def context_for_record(
    record: SiteRecord,
    annotation: AnnotationIndex,
) -> str:
    if record.context == "unannotated":
        record.context = annotation.classify(record.key)
    return record.context


def make_case_analysis(
    *,
    case: CaseSite,
    case_record: SiteRecord,
    case_score: SiteScore,
    local_records: Sequence[SiteRecord],
    global_records: Sequence[SiteRecord],
    score_cache: Mapping[tuple[str, SiteKey, int], Optional[SiteScore]],
    reference: SequenceReference,
    composition_radius: int,
    min_support: int,
    treatment_to_index: Mapping[str, int],
    min_local_controls: int,
    min_global_controls: int,
    domain_endpoint: str,
    matches_per_case: int,
    min_matches_per_case: int,
    pl_caliper: float,
    log_depth_caliper: float,
    log_support_caliper: float,
    gc_caliper: float,
    target_caliper: float,
) -> Optional[CaseAnalysis]:
    """
    Construct two independent case-specific matched sets:

    - local controls from the same genomic neighbourhood;
    - global controls from the random callable reservoir.

    Exact eligibility is handled before this function. Matching is based only
    on pre-treatment/technical covariates and never on entropy outcomes.
    """
    indices = required_indices(case, treatment_to_index)
    case_pl, case_depth = case_covariates(
        case,
        case_record,
        treatment_to_index,
    )
    case_gc, case_target_fraction = reference.local_composition(
        case.key,
        case_record.name,
        composition_radius,
    )
    if (
        case_pl is None
        or case_depth is None
        or not math.isfinite(case_score.site_support)
        or not math.isfinite(case_gc)
        or not math.isfinite(case_target_fraction)
    ):
        return None

    def prepare(
        records: Sequence[SiteRecord],
        endpoint: str,
        minimum_candidate_pool: int,
    ) -> Optional[
        tuple[
            np.ndarray,
            np.ndarray,
            tuple[SiteKey, ...],
            float,
            int,
            float,
            float,
            float,
            float,
            float,
        ]
    ]:
        values: list[float] = []
        keys: list[SiteKey] = []
        pls: list[float] = []
        depths: list[float] = []
        supports: list[float] = []
        gc_values: list[float] = []
        target_values: list[float] = []
        base_weights: list[float] = []

        for record in records:
            if (
                case_record.name in DNA_BASES
                and record.name in DNA_BASES
                and record.name != case_record.name
            ):
                continue

            control_pl = relevant_average(record.pls, indices)
            depth = relevant_average(record.depths, indices)
            score = score_cache.get(
                (case.mod_type, record.key, min_support)
            )
            if (
                control_pl is None
                or depth is None
                or score is None
            ):
                continue

            value = score.endpoint(endpoint)
            if not math.isfinite(value):
                continue

            gc_fraction, target_fraction = reference.local_composition(
                record.key,
                record.name,
                composition_radius,
            )
            if (
                not math.isfinite(gc_fraction)
                or not math.isfinite(target_fraction)
                or not math.isfinite(score.site_support)
            ):
                continue

            values.append(value)
            keys.append(record.key)
            pls.append(control_pl)
            depths.append(depth)
            supports.append(score.site_support)
            gc_values.append(gc_fraction)
            target_values.append(target_fraction)
            base_weights.append(record.sampling_weight)

        if len(values) < minimum_candidate_pool:
            return None

        values_array = np.asarray(values, dtype=float)
        pls_array = np.asarray(pls, dtype=float)
        depths_array = np.asarray(depths, dtype=float)
        supports_array = np.asarray(supports, dtype=float)
        gc_array = np.asarray(gc_values, dtype=float)
        target_array = np.asarray(target_values, dtype=float)
        base_array = np.asarray(base_weights, dtype=float)

        weights, match_mismatch = strict_match_weights(
            case_pl,
            case_depth,
            case_score.site_support,
            case_gc,
            case_target_fraction,
            pls_array,
            depths_array,
            supports_array,
            gc_array,
            target_array,
            base_array,
            matches_per_case=matches_per_case,
            min_matches_per_case=min_matches_per_case,
            pl_caliper=pl_caliper,
            log_depth_caliper=log_depth_caliper,
            log_support_caliper=log_support_caliper,
            gc_caliper=gc_caliper,
            target_caliper=target_caliper,
        )

        selected = np.isfinite(weights) & (weights > 0)
        if int(np.sum(selected)) < min_matches_per_case:
            return None

        values_array = values_array[selected]
        weights = weights[selected]
        selected_indices = np.flatnonzero(selected)
        selected_keys = tuple(keys[index] for index in selected_indices)
        pls_array = pls_array[selected]
        depths_array = depths_array[selected]
        supports_array = supports_array[selected]
        gc_array = gc_array[selected]
        target_array = target_array[selected]

        return (
            values_array,
            weights,
            selected_keys,
            match_mismatch,
            len(values_array),
            weighted_mean(pls_array, weights),
            weighted_mean(depths_array, weights),
            weighted_mean(supports_array, weights),
            weighted_mean(gc_array, weights),
            weighted_mean(target_array, weights),
        )

    local_prepared = prepare(
        local_records,
        "site_me",
        min_local_controls,
    )
    global_prepared = prepare(
        global_records,
        domain_endpoint,
        min_global_controls,
    )
    if local_prepared is None or global_prepared is None:
        return None

    (
        local_values,
        local_weights,
        local_control_keys,
        local_match_mismatch,
        local_n,
        local_expected_pl,
        local_expected_depth,
        local_expected_support,
        local_expected_gc,
        local_expected_target,
    ) = local_prepared
    (
        global_values,
        global_weights,
        global_control_keys,
        global_match_mismatch,
        global_n,
        global_expected_pl,
        global_expected_depth,
        global_expected_support,
        global_expected_gc,
        global_expected_target,
    ) = global_prepared

    # Mean matched-control contrasts are symmetric under within-set
    # randomization. This is preferable to mixing a case value with a median
    # reference when the formal null permutes labels within the matched set.
    local_reference = weighted_mean(
        local_values,
        local_weights,
    )
    global_reference = weighted_mean(
        global_values,
        global_weights,
    )

    local_percentile = (
        float(np.sum(
            local_weights[local_values <= case_score.site_me]
        ))
        / float(np.sum(local_weights))
    )
    case_domain = case_score.endpoint(domain_endpoint)
    global_percentile = (
        float(np.sum(
            global_weights[global_values <= case_domain]
        ))
        / float(np.sum(global_weights))
    )

    return CaseAnalysis(
        case=case,
        case_record=case_record,
        case_score=case_score,
        local_control_values=local_values,
        local_control_weights=local_weights,
        local_control_keys=local_control_keys,
        global_control_values=global_values,
        global_control_weights=global_weights,
        global_control_keys=global_control_keys,
        local_effect=case_score.site_me - local_reference,
        local_percentile=local_percentile,
        global_effect=case_domain - global_reference,
        global_percentile=global_percentile,
        local_ess=effective_sample_size(local_weights),
        global_ess=effective_sample_size(global_weights),
        local_controls=local_n,
        global_controls=global_n,
        # Maximum absolute standardized mean covariate mismatch within the
        # selected case-specific control set.
        local_within_case_max_standardized_mismatch=local_match_mismatch,
        global_within_case_max_standardized_mismatch=global_match_mismatch,
        case_control_pl=case_pl,
        case_control_depth=case_depth,
        case_me_support=case_score.site_support,
        case_gc_fraction=case_gc,
        case_target_base_fraction=case_target_fraction,
        local_expected_control_pl=local_expected_pl,
        local_expected_control_depth=local_expected_depth,
        local_expected_me_support=local_expected_support,
        local_expected_gc_fraction=local_expected_gc,
        local_expected_target_base_fraction=local_expected_target,
        global_expected_control_pl=global_expected_pl,
        global_expected_control_depth=global_expected_depth,
        global_expected_me_support=global_expected_support,
        global_expected_gc_fraction=global_expected_gc,
        global_expected_target_base_fraction=global_expected_target,
    )


def build_analyses_for_configuration(
    *,
    connection: sqlite3.Connection,
    mod_type: str,
    cases: Mapping[SiteKey, CaseSite],
    case_records: Mapping[tuple[str, SiteKey], SiteRecord],
    all_case_keys: set[SiteKey],
    reservoirs: Mapping[tuple[str, str, int], Reservoir],
    score_cache: dict[tuple[str, SiteKey, int], Optional[SiteScore]],
    me_index: MEIndex,
    annotation: AnnotationIndex,
    reference: SequenceReference,
    composition_radius: int,
    treatment_count: int,
    treatment_to_index: Mapping[str, int],
    min_support: int,
    local_radius: int,
    domain_radius: int,
    min_local_controls: int,
    min_global_controls: int,
    matches_per_case: int,
    min_matches_per_case: int,
    pl_caliper: float,
    log_depth_caliper: float,
    log_support_caliper: float,
    gc_caliper: float,
    target_caliper: float,
    context_match_global: bool = False,
) -> list[CaseAnalysis]:
    """
    Rebuild the exact case-specific eligibility graph for one sensitivity
    configuration. This avoids reusing a control merely because it was eligible
    for another case.
    """
    analyses: list[CaseAnalysis] = []
    domain_endpoint = f"domain_{domain_radius}"

    for key, case in cases.items():
        record = case_records.get((mod_type, key))
        if record is None:
            continue

        case_cache_key = (mod_type, key, min_support)
        if case_cache_key not in score_cache:
            score_cache[case_cache_key] = me_index.score(
                mod_type,
                key,
                min_support,
            )
        case_score = score_cache[case_cache_key]
        if case_score is None:
            continue

        local_records = query_local_controls(
            connection,
            mod_type,
            case,
            treatment_count,
            treatment_to_index,
            local_radius,
            all_case_keys,
        )
        for control in local_records:
            cache_key = (mod_type, control.key, min_support)
            if cache_key not in score_cache:
                score_cache[cache_key] = me_index.score(
                    mod_type,
                    control.key,
                    min_support,
                )

        global_records = eligible_global_controls(
            reservoirs,
            case,
            treatment_to_index,
        )
        if context_match_global:
            case_context = annotation.classify(case.key)
            global_records = [
                control
                for control in global_records
                if context_for_record(control, annotation) == case_context
            ]

        analysis = make_case_analysis(
            case=case,
            case_record=record,
            case_score=case_score,
            local_records=local_records,
            global_records=global_records,
            score_cache=score_cache,
            reference=reference,
            composition_radius=composition_radius,
            min_support=min_support,
            treatment_to_index=treatment_to_index,
            min_local_controls=min_local_controls,
            min_global_controls=min_global_controls,
            domain_endpoint=domain_endpoint,
            matches_per_case=matches_per_case,
            min_matches_per_case=min_matches_per_case,
            pl_caliper=pl_caliper,
            log_depth_caliper=log_depth_caliper,
            log_support_caliper=log_support_caliper,
            gc_caliper=gc_caliper,
            target_caliper=target_caliper,
        )
        if analysis is not None:
            analyses.append(analysis)

    return analyses


def standardized_mean_difference(
    case_values: np.ndarray,
    control_values: np.ndarray,
    scale_floor: float = 0.0,
) -> float:
    mask = np.isfinite(case_values) & np.isfinite(control_values)
    case_values = case_values[mask]
    control_values = control_values[mask]
    if len(case_values) < 2:
        return float("nan")
    pooled = math.sqrt(
        (
            float(np.var(case_values, ddof=1))
            + float(np.var(control_values, ddof=1))
        )
        / 2.0
    )
    denominator = max(pooled, scale_floor)
    if denominator <= 0:
        return 0.0
    return (
        float(np.mean(case_values))
        - float(np.mean(control_values))
    ) / denominator


def make_balance_table(
    analyses_by_mod: Mapping[str, Sequence[CaseAnalysis]],
) -> pd.DataFrame:
    rows = []
    for mod_type, analyses in analyses_by_mod.items():
        if not analyses:
            continue
        case_arrays = {
            "control_pl": np.asarray(
                [item.case_control_pl for item in analyses],
                dtype=float,
            ),
            "log_control_depth": np.log1p(
                np.asarray(
                    [item.case_control_depth for item in analyses],
                    dtype=float,
                )
            ),
            "log_me_support": np.log1p(
                np.asarray(
                    [item.case_me_support for item in analyses],
                    dtype=float,
                )
            ),
            "gc_fraction": np.asarray(
                [item.case_gc_fraction for item in analyses],
                dtype=float,
            ),
            "target_base_fraction": np.asarray(
                [
                    item.case_target_base_fraction
                    for item in analyses
                ],
                dtype=float,
            ),
        }
        for endpoint in ("local", "global"):
            control_arrays = {
                "control_pl": np.asarray(
                    [
                        (
                            item.local_expected_control_pl
                            if endpoint == "local"
                            else item.global_expected_control_pl
                        )
                        for item in analyses
                    ],
                    dtype=float,
                ),
                "log_control_depth": np.log1p(
                    np.asarray(
                        [
                            (
                                item.local_expected_control_depth
                                if endpoint == "local"
                                else item.global_expected_control_depth
                            )
                            for item in analyses
                        ],
                        dtype=float,
                    )
                ),
                "log_me_support": np.log1p(
                    np.asarray(
                        [
                            (
                                item.local_expected_me_support
                                if endpoint == "local"
                                else item.global_expected_me_support
                            )
                            for item in analyses
                        ],
                        dtype=float,
                    )
                ),
                "gc_fraction": np.asarray(
                    [
                        (
                            item.local_expected_gc_fraction
                            if endpoint == "local"
                            else item.global_expected_gc_fraction
                        )
                        for item in analyses
                    ],
                    dtype=float,
                ),
                "target_base_fraction": np.asarray(
                    [
                        (
                            item.local_expected_target_base_fraction
                            if endpoint == "local"
                            else item.global_expected_target_base_fraction
                        )
                        for item in analyses
                    ],
                    dtype=float,
                ),
            }
            scale_floors = {
                "control_pl": 0.05,
                "log_control_depth": 0.25,
                "log_me_support": 0.25,
                "gc_fraction": 0.02,
                "target_base_fraction": 0.02,
            }
            for covariate in case_arrays:
                case_values = case_arrays[covariate]
                control_values = control_arrays[covariate]
                rows.append(
                    {
                        "mod_type": mod_type,
                        "endpoint": endpoint,
                        "covariate": covariate,
                        "case_mean": float(np.nanmean(case_values)),
                        "matched_control_mean": float(
                            np.nanmean(control_values)
                        ),
                        "standardized_mean_difference": (
                            standardized_mean_difference(
                                case_values,
                                control_values,
                                scale_floor=scale_floors[covariate],
                            )
                        ),
                        "mean_absolute_pair_difference": float(
                            np.nanmean(
                                np.abs(case_values - control_values)
                            )
                        ),
                        "standardization_scale_floor": (
                            scale_floors[covariate]
                        ),
                        "case_count": len(analyses),
                    }
                )
    return pd.DataFrame(rows)


# =============================================================================
# RANDOMIZATION AND BOOTSTRAP
# =============================================================================

def spatial_clusters(
    analyses: Sequence[CaseAnalysis],
    cluster_bp: int,
) -> list[list[int]]:
    order = sorted(
        range(len(analyses)),
        key=lambda index: (
            analyses[index].case.key.chrom,
            analyses[index].case.key.start,
            analyses[index].case.key.strand,
        ),
    )
    clusters: list[list[int]] = []
    for index in order:
        analysis = analyses[index]
        if not clusters:
            clusters.append([index])
            continue
        previous = analyses[clusters[-1][-1]]
        if (
            analysis.case.key.chrom == previous.case.key.chrom
            and analysis.case.key.strand == previous.case.key.strand
            and (
                analysis.case.key.start
                - previous.case.key.start
                <= cluster_bp
            )
        ):
            clusters[-1].append(index)
        else:
            clusters.append([index])
    return clusters


def aggregate_by_clusters(
    values: np.ndarray,
    clusters: Sequence[Sequence[int]],
) -> np.ndarray:
    return np.asarray(
        [
            float(np.mean(values[list(cluster)]))
            for cluster in clusters
        ],
        dtype=float,
    )


def matched_covariate_arrays(
    analyses: Sequence[CaseAnalysis],
    endpoint: str,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    case_arrays = {
        "control_pl": np.asarray(
            [item.case_control_pl for item in analyses],
            dtype=float,
        ),
        "log_control_depth": np.log1p(np.asarray(
            [item.case_control_depth for item in analyses],
            dtype=float,
        )),
        "log_me_support": np.log1p(np.asarray(
            [item.case_me_support for item in analyses],
            dtype=float,
        )),
        "gc_fraction": np.asarray(
            [item.case_gc_fraction for item in analyses],
            dtype=float,
        ),
        "target_base_fraction": np.asarray(
            [item.case_target_base_fraction for item in analyses],
            dtype=float,
        ),
    }
    if endpoint == "local":
        control_arrays = {
            "control_pl": np.asarray(
                [item.local_expected_control_pl for item in analyses],
                dtype=float,
            ),
            "log_control_depth": np.log1p(np.asarray(
                [item.local_expected_control_depth for item in analyses],
                dtype=float,
            )),
            "log_me_support": np.log1p(np.asarray(
                [item.local_expected_me_support for item in analyses],
                dtype=float,
            )),
            "gc_fraction": np.asarray(
                [item.local_expected_gc_fraction for item in analyses],
                dtype=float,
            ),
            "target_base_fraction": np.asarray(
                [
                    item.local_expected_target_base_fraction
                    for item in analyses
                ],
                dtype=float,
            ),
        }
    else:
        control_arrays = {
            "control_pl": np.asarray(
                [item.global_expected_control_pl for item in analyses],
                dtype=float,
            ),
            "log_control_depth": np.log1p(np.asarray(
                [item.global_expected_control_depth for item in analyses],
                dtype=float,
            )),
            "log_me_support": np.log1p(np.asarray(
                [item.global_expected_me_support for item in analyses],
                dtype=float,
            )),
            "gc_fraction": np.asarray(
                [item.global_expected_gc_fraction for item in analyses],
                dtype=float,
            ),
            "target_base_fraction": np.asarray(
                [
                    item.global_expected_target_base_fraction
                    for item in analyses
                ],
                dtype=float,
            ),
        }
    return case_arrays, control_arrays


BALANCE_SCALE_FLOORS = {
    "control_pl": 0.05,
    "log_control_depth": 0.25,
    "log_me_support": 0.25,
    "gc_fraction": 0.02,
    "target_base_fraction": 0.02,
}


def matched_balance_metrics(
    analyses: Sequence[CaseAnalysis],
    endpoint: str,
) -> tuple[dict[str, float], float]:
    case_arrays, control_arrays = matched_covariate_arrays(
        analyses,
        endpoint,
    )
    smds = {
        f"smd_{name}": standardized_mean_difference(
            case_arrays[name],
            control_arrays[name],
            scale_floor=BALANCE_SCALE_FLOORS[name],
        )
        for name in case_arrays
    }
    finite = [
        abs(value)
        for value in smds.values()
        if math.isfinite(value)
    ]
    return smds, (max(finite) if finite else float("nan"))


def trim_for_covariate_balance(
    analyses: Sequence[CaseAnalysis],
    endpoint: str,
    *,
    threshold: float = 0.10,
    minimum_cases: int = 10,
    maximum_trim_fraction: float = 0.50,
) -> tuple[list[CaseAnalysis], dict[str, object]]:
    """
    Outcome-blind common-support trimming.

    Strict calipers are applied first. If small directional residual imbalance
    remains across the matched sets, cases contributing most to the worst
    covariate imbalance are removed without consulting ME outcomes. Trimming
    stops at the prespecified SMD threshold or the maximum allowed attrition.
    """
    retained = list(analyses)
    removed_case_ids: list[str] = []
    original_count = len(retained)
    minimum_retained = max(
        minimum_cases,
        int(math.ceil(
            original_count * (1.0 - maximum_trim_fraction)
        )),
    )

    while len(retained) > minimum_retained:
        smds, maximum = matched_balance_metrics(retained, endpoint)
        if math.isfinite(maximum) and maximum <= threshold:
            break
        if not math.isfinite(maximum):
            break

        worst_key = max(
            smds,
            key=lambda key: (
                abs(smds[key])
                if math.isfinite(smds[key])
                else -1.0
            ),
        )
        covariate = worst_key.replace("smd_", "", 1)
        case_arrays, control_arrays = matched_covariate_arrays(
            retained,
            endpoint,
        )
        case_values = case_arrays[covariate]
        control_values = control_arrays[covariate]
        pooled = math.sqrt(
            (
                float(np.var(case_values, ddof=1))
                + float(np.var(control_values, ddof=1))
            )
            / 2.0
        )
        scale = max(
            pooled,
            BALANCE_SCALE_FLOORS[covariate],
        )
        signed_pair_difference = (
            case_values - control_values
        ) / scale
        direction = 1.0 if smds[worst_key] >= 0 else -1.0
        contribution = direction * signed_pair_difference
        contribution[~np.isfinite(contribution)] = np.inf
        remove_index = int(np.argmax(contribution))
        removed = retained.pop(remove_index)
        removed_case_ids.append(
            f"{removed.case.mod_type}|{removed.case.key.chrom}:"
            f"{removed.case.key.start}-{removed.case.key.end}"
            f"({removed.case.key.strand})"
        )

    final_smds, final_maximum = matched_balance_metrics(
        retained,
        endpoint,
    )
    return retained, {
        "site_cases_before_balance_trim": original_count,
        "balance_trimmed_cases": original_count - len(retained),
        "balance_trimmed_case_ids": ",".join(removed_case_ids),
        "balance_trim_fraction": (
            (original_count - len(retained)) / original_count
            if original_count
            else float("nan")
        ),
        "balance_trim_limit_fraction": maximum_trim_fraction,
        "balance_target_max_abs_smd": threshold,
        "balance_target_achieved": (
            math.isfinite(final_maximum)
            and final_maximum <= threshold
        ),
        **final_smds,
        "max_abs_covariate_smd": final_maximum,
    }


def randomization_test(
    analyses: Sequence[CaseAnalysis],
    endpoint: str,
    *,
    randomizations: int,
    bootstraps: int,
    cluster_bp: int,
    seed: int,
) -> dict[str, object]:
    """
    Exact conditional label randomization within case-specific matched sets.

    Each set contains one observed case and its strictly calipered controls.
    Under the matched-set null, any member could have carried the case label.
    Spatially adjacent cases are first collapsed to event-cluster means before
    the genome-level statistic is calculated.
    """
    if endpoint not in {"local", "global"}:
        raise ValueError(endpoint)

    analyses, balance_trim = trim_for_covariate_balance(
        analyses,
        endpoint,
        threshold=0.10,
        minimum_cases=10,
        maximum_trim_fraction=0.50,
    )

    observed_site_effects = np.asarray(
        [
            analysis.local_effect
            if endpoint == "local"
            else analysis.global_effect
            for analysis in analyses
        ],
        dtype=float,
    )
    clusters = spatial_clusters(analyses, cluster_bp)
    observed_cluster = aggregate_by_clusters(
        observed_site_effects,
        clusters,
    )
    observed = (
        float(np.mean(observed_cluster))
        if len(observed_cluster)
        else float("nan")
    )

    rng = np.random.default_rng(seed)
    null = np.zeros(randomizations, dtype=float)
    null_tail90 = np.zeros(randomizations, dtype=float)
    null_tail95 = np.zeros(randomizations, dtype=float)

    for cluster in clusters:
        cluster_effect = np.zeros(randomizations, dtype=float)
        cluster_tail90 = np.zeros(randomizations, dtype=float)
        cluster_tail95 = np.zeros(randomizations, dtype=float)

        for index in cluster:
            analysis = analyses[index]
            if endpoint == "local":
                case_value = float(analysis.case_score.site_me)
                controls = np.asarray(
                    analysis.local_control_values,
                    dtype=float,
                )
            else:
                domain_name = None
                # The stored global effect and matched control mean recover the
                # exact case-domain value without relying on one fixed radius.
                controls = np.asarray(
                    analysis.global_control_values,
                    dtype=float,
                )
                case_value = float(
                    analysis.global_effect
                    + np.mean(controls)
                )

            members = np.concatenate(
                [np.asarray([case_value], dtype=float), controls]
            )
            member_count = len(members)
            if member_count < 2:
                continue

            chosen = rng.integers(
                0,
                member_count,
                size=randomizations,
            )
            chosen_values = members[chosen]
            other_means = (
                np.sum(members) - chosen_values
            ) / float(member_count - 1)
            cluster_effect += chosen_values - other_means

            pseudo_percentiles = np.empty(member_count, dtype=float)
            for member_index, member_value in enumerate(members):
                other_mask = np.ones(member_count, dtype=bool)
                other_mask[member_index] = False
                others = members[other_mask]
                pseudo_percentiles[member_index] = float(
                    np.mean(others <= member_value)
                )
            sampled_percentiles = pseudo_percentiles[chosen]
            cluster_tail90 += sampled_percentiles >= 0.90
            cluster_tail95 += sampled_percentiles >= 0.95

        cluster_size = float(len(cluster))
        null += cluster_effect / cluster_size
        null_tail90 += cluster_tail90 / cluster_size
        null_tail95 += cluster_tail95 / cluster_size

    cluster_count = float(len(clusters))
    if cluster_count > 0:
        null /= cluster_count
        null_tail90 /= cluster_count
        null_tail95 /= cluster_count

    p_high = (
        int(np.sum(null >= observed)) + 1
    ) / (randomizations + 1)
    p_low = (
        int(np.sum(null <= observed)) + 1
    ) / (randomizations + 1)
    p_two = (
        int(np.sum(np.abs(null) >= abs(observed))) + 1
    ) / (randomizations + 1)

    bootstrap = np.empty(bootstraps, dtype=float)
    if len(observed_cluster):
        for iteration in range(bootstraps):
            sampled = rng.integers(
                0,
                len(observed_cluster),
                size=len(observed_cluster),
            )
            bootstrap[iteration] = float(
                np.mean(observed_cluster[sampled])
            )
    else:
        bootstrap.fill(np.nan)

    percentiles = np.asarray(
        [
            analysis.local_percentile
            if endpoint == "local"
            else analysis.global_percentile
            for analysis in analyses
        ],
        dtype=float,
    )
    observed_tail90 = (
        float(np.mean(
            aggregate_by_clusters(
                (percentiles >= 0.90).astype(float),
                clusters,
            )
        ))
        if clusters
        else float("nan")
    )
    observed_tail95 = (
        float(np.mean(
            aggregate_by_clusters(
                (percentiles >= 0.95).astype(float),
                clusters,
            )
        ))
        if clusters
        else float("nan")
    )

    case_pl = np.asarray(
        [item.case_control_pl for item in analyses],
        dtype=float,
    )
    case_depth = np.log1p(
        np.asarray(
            [item.case_control_depth for item in analyses],
            dtype=float,
        )
    )
    case_support = np.log1p(
        np.asarray(
            [item.case_me_support for item in analyses],
            dtype=float,
        )
    )
    case_gc = np.asarray(
        [item.case_gc_fraction for item in analyses],
        dtype=float,
    )
    case_target = np.asarray(
        [item.case_target_base_fraction for item in analyses],
        dtype=float,
    )

    if endpoint == "local":
        control_pl = np.asarray(
            [item.local_expected_control_pl for item in analyses],
            dtype=float,
        )
        control_depth = np.log1p(
            np.asarray(
                [
                    item.local_expected_control_depth
                    for item in analyses
                ],
                dtype=float,
            )
        )
        control_support = np.log1p(
            np.asarray(
                [
                    item.local_expected_me_support
                    for item in analyses
                ],
                dtype=float,
            )
        )
        control_gc = np.asarray(
            [
                item.local_expected_gc_fraction
                for item in analyses
            ],
            dtype=float,
        )
        control_target = np.asarray(
            [
                item.local_expected_target_base_fraction
                for item in analyses
            ],
            dtype=float,
        )
        matched_counts = np.asarray(
            [item.local_controls for item in analyses],
            dtype=float,
        )
        within_case_mismatch = np.asarray(
            [
                item.local_within_case_max_standardized_mismatch
                for item in analyses
            ],
            dtype=float,
        )
    else:
        control_pl = np.asarray(
            [item.global_expected_control_pl for item in analyses],
            dtype=float,
        )
        control_depth = np.log1p(
            np.asarray(
                [
                    item.global_expected_control_depth
                    for item in analyses
                ],
                dtype=float,
            )
        )
        control_support = np.log1p(
            np.asarray(
                [
                    item.global_expected_me_support
                    for item in analyses
                ],
                dtype=float,
            )
        )
        control_gc = np.asarray(
            [
                item.global_expected_gc_fraction
                for item in analyses
            ],
            dtype=float,
        )
        control_target = np.asarray(
            [
                item.global_expected_target_base_fraction
                for item in analyses
            ],
            dtype=float,
        )
        matched_counts = np.asarray(
            [item.global_controls for item in analyses],
            dtype=float,
        )
        within_case_mismatch = np.asarray(
            [
                item.global_within_case_max_standardized_mismatch
                for item in analyses
            ],
            dtype=float,
        )

    balance_smds = {
        "smd_control_pl": standardized_mean_difference(
            case_pl,
            control_pl,
            scale_floor=0.05,
        ),
        "smd_log_control_depth": standardized_mean_difference(
            case_depth,
            control_depth,
            scale_floor=0.25,
        ),
        "smd_log_me_support": standardized_mean_difference(
            case_support,
            control_support,
            scale_floor=0.25,
        ),
        "smd_gc_fraction": standardized_mean_difference(
            case_gc,
            control_gc,
            scale_floor=0.02,
        ),
        "smd_target_base_fraction": standardized_mean_difference(
            case_target,
            control_target,
            scale_floor=0.02,
        ),
    }
    finite_balance = [
        abs(value)
        for value in balance_smds.values()
        if math.isfinite(value)
    ]
    max_abs_smd = (
        max(finite_balance)
        if finite_balance
        else float("nan")
    )

    if len(analyses) < 10:
        validity = "underpowered_fewer_than_10_valid_cases"
    elif len(clusters) < 5:
        validity = "underpowered_too_few_spatial_clusters"
    elif not math.isfinite(max_abs_smd):
        validity = "invalid_missing_balance_diagnostics"
    elif max_abs_smd > 0.10:
        validity = "invalid_poor_covariate_balance"
    else:
        validity = "valid_matched_set_genomic_test"

    return {
        "endpoint": endpoint,
        "site_cases": len(analyses),
        "event_clusters": len(clusters),
        "cluster_bp": cluster_bp,
        "minimum_cluster_rule_passed": len(clusters) >= 5,
        "inference_validity": validity,
        **balance_trim,
        "max_abs_covariate_smd": max_abs_smd,
        **balance_smds,
        "median_controls_per_case": (
            float(np.median(matched_counts))
            if len(matched_counts)
            else float("nan")
        ),
        "minimum_controls_per_case": (
            int(np.min(matched_counts))
            if len(matched_counts)
            else 0
        ),
        "median_within_case_standardized_mismatch": (
            float(np.nanmedian(within_case_mismatch))
            if len(within_case_mismatch)
            else float("nan")
        ),
        "observed_mean_effect": observed,
        "bootstrap_ci_low": (
            float(np.nanquantile(bootstrap, 0.025))
            if np.any(np.isfinite(bootstrap))
            else float("nan")
        ),
        "bootstrap_ci_high": (
            float(np.nanquantile(bootstrap, 0.975))
            if np.any(np.isfinite(bootstrap))
            else float("nan")
        ),
        "genomic_randomization_p_high": p_high,
        "genomic_randomization_p_low": p_low,
        "genomic_randomization_p_two_sided": p_two,
        "null_mean": float(np.mean(null)),
        "null_sd": float(np.std(null, ddof=1)),
        "median_case_percentile": (
            float(np.median(percentiles))
            if len(percentiles)
            else float("nan")
        ),
        "fraction_case_percentile_ge_0.90": observed_tail90,
        "fraction_case_percentile_ge_0.95": observed_tail95,
        "tail90_randomization_p_high": (
            (int(np.sum(null_tail90 >= observed_tail90)) + 1)
            / (randomizations + 1)
            if math.isfinite(observed_tail90)
            else float("nan")
        ),
        "tail95_randomization_p_high": (
            (int(np.sum(null_tail95 >= observed_tail95)) + 1)
            / (randomizations + 1)
            if math.isfinite(observed_tail95)
            else float("nan")
        ),
        "randomizations": randomizations,
        "bootstraps": bootstraps,
        "inference_scope": (
            "case_specific_caliper_matched_set_genomic_randomization"
        ),
    }


def direction_difference_test(
    gain: Sequence[CaseAnalysis],
    loss: Sequence[CaseAnalysis],
    endpoint: str,
    *,
    randomizations: int,
    bootstraps: int,
    cluster_bp: int,
    seed: int,
) -> dict[str, object]:
    """
    Exploratory comparison of matched-set ME effects for PL gains versus losses.

    Each direction is first trimmed independently for covariate balance. The
    null distribution is generated by exact within-set label randomization
    inside each direction, then differenced. This is not a randomized
    gain/loss assignment and should remain exploratory.
    """
    gain, gain_trim = trim_for_covariate_balance(
        gain,
        endpoint,
        threshold=0.10,
        minimum_cases=5,
        maximum_trim_fraction=0.50,
    )
    loss, loss_trim = trim_for_covariate_balance(
        loss,
        endpoint,
        threshold=0.10,
        minimum_cases=5,
        maximum_trim_fraction=0.50,
    )
    if len(gain) < 5 or len(loss) < 5:
        return {
            "endpoint": endpoint,
            "gain_cases": len(gain),
            "loss_cases": len(loss),
            "status": "untestable_insufficient_group",
            "gain_balance_trimmed_cases": gain_trim.get(
                "balance_trimmed_cases", 0
            ),
            "loss_balance_trimmed_cases": loss_trim.get(
                "balance_trimmed_cases", 0
            ),
        }

    def effects(group: Sequence[CaseAnalysis]) -> np.ndarray:
        return np.asarray(
            [
                item.local_effect
                if endpoint == "local"
                else item.global_effect
                for item in group
            ],
            dtype=float,
        )

    def group_null(
        group: Sequence[CaseAnalysis],
        clusters: Sequence[Sequence[int]],
        rng: np.random.Generator,
    ) -> np.ndarray:
        null = np.zeros(randomizations, dtype=float)
        for cluster in clusters:
            cluster_effect = np.zeros(randomizations, dtype=float)
            for index in cluster:
                analysis = group[index]
                if endpoint == "local":
                    case_value = float(analysis.case_score.site_me)
                    controls = np.asarray(
                        analysis.local_control_values,
                        dtype=float,
                    )
                else:
                    controls = np.asarray(
                        analysis.global_control_values,
                        dtype=float,
                    )
                    case_value = float(
                        analysis.global_effect
                        + np.mean(controls)
                    )
                members = np.concatenate(
                    [np.asarray([case_value]), controls]
                )
                chosen = rng.integers(
                    0,
                    len(members),
                    size=randomizations,
                )
                chosen_values = members[chosen]
                other_means = (
                    np.sum(members) - chosen_values
                ) / float(len(members) - 1)
                cluster_effect += chosen_values - other_means
            null += cluster_effect / float(len(cluster))
        return null / float(len(clusters))

    gain_clusters = spatial_clusters(gain, cluster_bp)
    loss_clusters = spatial_clusters(loss, cluster_bp)
    gain_effect = aggregate_by_clusters(
        effects(gain),
        gain_clusters,
    )
    loss_effect = aggregate_by_clusters(
        effects(loss),
        loss_clusters,
    )
    observed = float(
        np.mean(gain_effect) - np.mean(loss_effect)
    )

    rng = np.random.default_rng(seed)
    null = (
        group_null(gain, gain_clusters, rng)
        - group_null(loss, loss_clusters, rng)
    )
    p_two = (
        int(np.sum(np.abs(null) >= abs(observed))) + 1
    ) / (randomizations + 1)

    bootstrap = np.empty(bootstraps, dtype=float)
    for iteration in range(bootstraps):
        gain_sample = gain_effect[
            rng.integers(
                0,
                len(gain_effect),
                size=len(gain_effect),
            )
        ]
        loss_sample = loss_effect[
            rng.integers(
                0,
                len(loss_effect),
                size=len(loss_effect),
            )
        ]
        bootstrap[iteration] = (
            float(np.mean(gain_sample))
            - float(np.mean(loss_sample))
        )

    return {
        "endpoint": endpoint,
        "gain_cases": len(gain),
        "loss_cases": len(loss),
        "gain_clusters": len(gain_clusters),
        "loss_clusters": len(loss_clusters),
        "gain_balance_trimmed_cases": gain_trim.get(
            "balance_trimmed_cases", 0
        ),
        "loss_balance_trimmed_cases": loss_trim.get(
            "balance_trimmed_cases", 0
        ),
        "gain_max_abs_covariate_smd": gain_trim.get(
            "max_abs_covariate_smd", float("nan")
        ),
        "loss_max_abs_covariate_smd": loss_trim.get(
            "max_abs_covariate_smd", float("nan")
        ),
        "gain_minus_loss_effect": observed,
        "bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
        "bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
        "randomization_p_two_sided": p_two,
        "status": "tested_exploratory_matched_within_direction",
        "inference_scope": (
            "exploratory_difference_of_within_direction_matched_set_effects"
        ),
    }


# =============================================================================
# CONTINUOUS RESPONSE MODELS
# =============================================================================

def two_sided_t_pvalue(statistic: float, degrees_of_freedom: int) -> float:
    if not math.isfinite(statistic) or degrees_of_freedom < 1:
        return float("nan")
    try:
        from scipy.stats import t as student_t  # type: ignore
        return float(
            2.0 * student_t.sf(abs(statistic), degrees_of_freedom)
        )
    except Exception:
        # Conservative normal approximation when SciPy is unavailable.
        return float(math.erfc(abs(statistic) / math.sqrt(2.0)))


def t_critical_975(degrees_of_freedom: int) -> float:
    if degrees_of_freedom < 1:
        return 1.96
    try:
        from scipy.stats import t as student_t  # type: ignore
        return float(student_t.ppf(0.975, degrees_of_freedom))
    except Exception:
        return 1.96


def weighted_center_scale(
    values: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    mean = weighted_mean(values, weights)
    variance = weighted_mean((values - mean) ** 2, weights)
    scale = math.sqrt(max(variance, 1e-12))
    return (values - mean) / scale


def weighted_cluster_regression(
    frame: pd.DataFrame,
    *,
    predictor: str,
    include_treatment_fixed_effects: bool,
) -> dict[str, object]:
    """
    Weighted least squares with a CR1 genomic-block sandwich covariance.

    The outcome is absolute treatment-minus-control PL. Baseline PL is modeled
    flexibly with a cubic polynomial; depth, ME support, sequence composition,
    overlap-window count and callability are adjusted as nuisance covariates.
    The reported ME coefficient is the adjusted change in |delta PL| per 0.10
    increase in baseline ME.
    """
    required = [
        "abs_delta_pl",
        predictor,
        "control_pl",
        "control_depth",
        "me_support",
        "gc_fraction",
        "target_base_fraction",
        "overlap_windows",
        "callable_treatment_count",
        "sampling_weight",
        "genomic_block",
        "treatment",
    ]
    subset = frame.dropna(subset=required).copy()
    for column in required:
        if column in {"genomic_block", "treatment"}:
            continue
        subset = subset[np.isfinite(
            pd.to_numeric(subset[column], errors="coerce")
        )]
    subset = subset[
        pd.to_numeric(
            subset["sampling_weight"],
            errors="coerce",
        ) > 0
    ]

    if len(subset) < 200:
        return {
            "status": "underpowered_fewer_than_200_rows",
            "rows": len(subset),
        }

    y = subset["abs_delta_pl"].to_numpy(dtype=float)
    me = subset[predictor].to_numpy(dtype=float)
    weights = subset["sampling_weight"].to_numpy(dtype=float)
    weights = weights / float(np.mean(weights))

    pl = subset["control_pl"].to_numpy(dtype=float) - 0.5
    log_depth = np.log1p(
        subset["control_depth"].to_numpy(dtype=float)
    )
    log_support = np.log1p(
        subset["me_support"].to_numpy(dtype=float)
    )
    gc = subset["gc_fraction"].to_numpy(dtype=float)
    target = subset["target_base_fraction"].to_numpy(dtype=float)
    overlap = np.log1p(
        subset["overlap_windows"].to_numpy(dtype=float)
    )
    callable_count = subset[
        "callable_treatment_count"
    ].to_numpy(dtype=float)

    columns = [
        np.ones(len(subset), dtype=float),
        me,
        weighted_center_scale(pl, weights),
        weighted_center_scale(pl ** 2, weights),
        weighted_center_scale(pl ** 3, weights),
        weighted_center_scale(log_depth, weights),
        weighted_center_scale(log_support, weights),
        weighted_center_scale(gc, weights),
        weighted_center_scale(target, weights),
        weighted_center_scale(overlap, weights),
        weighted_center_scale(callable_count, weights),
    ]
    names = [
        "intercept",
        "baseline_ME",
        "baseline_PL_linear",
        "baseline_PL_quadratic",
        "baseline_PL_cubic",
        "log_control_depth",
        "log_ME_support",
        "GC_fraction",
        "target_base_fraction",
        "log_overlap_windows",
        "callable_treatment_count",
    ]

    if include_treatment_fixed_effects:
        treatments = sorted(
            str(value)
            for value in subset["treatment"].unique()
        )
        for treatment in treatments[1:]:
            columns.append(
                (
                    subset["treatment"].astype(str).to_numpy()
                    == treatment
                ).astype(float)
            )
            names.append(f"treatment[{treatment}]")

    X = np.column_stack(columns)
    n, p = X.shape
    weighted_X = X * weights[:, None]
    xtwx = X.T @ weighted_X
    xtwy = X.T @ (weights * y)
    bread = np.linalg.pinv(xtwx, rcond=1e-10)
    beta = bread @ xtwy
    residual = y - X @ beta

    cluster_labels = subset["genomic_block"].astype(str).to_numpy()
    unique_clusters = np.unique(cluster_labels)
    cluster_count = len(unique_clusters)
    meat = np.zeros((p, p), dtype=float)
    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        score = X[mask].T @ (
            weights[mask] * residual[mask]
        )
        meat += np.outer(score, score)

    covariance = bread @ meat @ bread
    if cluster_count > 1 and n > p:
        covariance *= (
            cluster_count / (cluster_count - 1.0)
        ) * (
            (n - 1.0) / (n - p)
        )

    coefficient = float(beta[1])
    variance = float(covariance[1, 1])
    standard_error = (
        math.sqrt(variance)
        if variance >= 0 and math.isfinite(variance)
        else float("nan")
    )
    statistic = (
        coefficient / standard_error
        if standard_error > 0
        else float("nan")
    )
    degrees_of_freedom = max(cluster_count - 1, 1)
    p_value = two_sided_t_pvalue(
        statistic,
        degrees_of_freedom,
    )
    critical = t_critical_975(degrees_of_freedom)

    weighted_y_mean = weighted_mean(y, weights)
    rss = float(np.sum(weights * residual ** 2))
    tss = float(np.sum(weights * (y - weighted_y_mean) ** 2))
    r_squared = (
        1.0 - rss / tss
        if tss > 0
        else float("nan")
    )
    adjusted_r_squared = (
        1.0 - (1.0 - r_squared) * (n - 1.0) / (n - p)
        if math.isfinite(r_squared) and n > p
        else float("nan")
    )
    partial_r_squared = (
        statistic ** 2
        / (statistic ** 2 + degrees_of_freedom)
        if math.isfinite(statistic)
        else float("nan")
    )

    weight_ess = effective_sample_size(weights)
    status = (
        "valid_cluster_robust_continuous_model"
        if (
            n >= 200
            and cluster_count >= 20
            and math.isfinite(standard_error)
            and standard_error > 0
        )
        else "underpowered_or_singular_continuous_model"
    )

    return {
        "status": status,
        "rows": n,
        "weighted_effective_rows": weight_ess,
        "genomic_blocks": cluster_count,
        "parameters": p,
        "predictor": predictor,
        "coefficient_per_1_ME": coefficient,
        "coefficient_per_0.1_ME": coefficient * 0.1,
        "percentage_point_abs_deltaPL_per_0.1_ME": (
            coefficient * 0.1 * 100.0
        ),
        "standard_error_per_1_ME": standard_error,
        "ci_low_per_0.1_ME": (
            (coefficient - critical * standard_error) * 0.1
            if math.isfinite(standard_error)
            else float("nan")
        ),
        "ci_high_per_0.1_ME": (
            (coefficient + critical * standard_error) * 0.1
            if math.isfinite(standard_error)
            else float("nan")
        ),
        "cluster_robust_t": statistic,
        "cluster_robust_df": degrees_of_freedom,
        "cluster_robust_p_two_sided": p_value,
        "weighted_R2": r_squared,
        "weighted_adjusted_R2": adjusted_r_squared,
        "partial_R2_baseline_ME": partial_r_squared,
        "mean_abs_delta_pl": weighted_y_mean,
        "mean_baseline_ME": weighted_mean(me, weights),
        "case_rows": int(subset["is_case"].astype(bool).sum()),
        "significant_case_rows": int(
            subset["significant_in_treatment"].astype(bool).sum()
        ),
        "model_terms": ",".join(names),
        "inference_scope": (
            "genomic_association_with_cluster_robust_blocks_not_"
            "biological_replicate_inference"
        ),
    }


def build_continuous_response_table(
    *,
    records_by_mod: Mapping[
        str,
        Mapping[SiteKey, SiteRecord],
    ],
    cases_by_mod: Mapping[
        str,
        Mapping[SiteKey, CaseSite],
    ],
    score_cache: Mapping[
        tuple[str, SiteKey, int],
        Optional[SiteScore],
    ],
    reference: SequenceReference,
    annotation: AnnotationIndex,
    treatments: Sequence[str],
    treatment_to_index: Mapping[str, int],
    min_support: int,
    domain_radius: int,
    composition_radius: int,
    block_bp: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    domain_endpoint = f"domain_{domain_radius}"

    for mod_type, record_map in records_by_mod.items():
        case_map = cases_by_mod.get(mod_type, {})
        for key, record in record_map.items():
            score = score_cache.get((mod_type, key, min_support))
            if score is None:
                continue
            domain_me = score.endpoint(domain_endpoint)
            if (
                not math.isfinite(score.site_me)
                or not math.isfinite(domain_me)
                or not math.isfinite(score.site_support)
            ):
                continue

            gc_fraction, target_fraction = reference.local_composition(
                key,
                record.name,
                composition_radius,
            )
            if (
                not math.isfinite(gc_fraction)
                or not math.isfinite(target_fraction)
            ):
                continue

            case = case_map.get(key)
            callable_count = bin(int(record.mask)).count("1")
            for treatment in treatments:
                index = treatment_to_index[treatment]
                control_pl = record.pls[index]
                treatment_pl = record.post_pls[index]
                control_depth = record.depths[index]
                if (
                    control_pl is None
                    or treatment_pl is None
                    or control_depth is None
                ):
                    continue
                delta = treatment_pl - control_pl
                rows.append(
                    {
                        "mod_type": mod_type,
                        "treatment": treatment,
                        "chrom": key.chrom,
                        "start": key.start,
                        "end": key.end,
                        "strand": key.strand,
                        "context": annotation.classify(key),
                        "genomic_block": (
                            f"{key.chrom}:{key.start // block_bp}"
                        ),
                        "control_pl": control_pl,
                        "treatment_pl": treatment_pl,
                        "delta_pl": delta,
                        "abs_delta_pl": abs(delta),
                        "control_depth": control_depth,
                        "site_me": score.site_me,
                        "domain_me": domain_me,
                        "context_500_me": (
                            score.context_500_excluding_site
                        ),
                        "local_prominence_500": (
                            score.local_prominence_500
                        ),
                        "local_percentile_500": (
                            score.local_percentile_500
                        ),
                        "local_peak90": (
                            float(score.local_percentile_500 >= 0.90)
                            if math.isfinite(score.local_percentile_500)
                            else float("nan")
                        ),
                        "local_peak95": (
                            float(score.local_percentile_500 >= 0.95)
                            if math.isfinite(score.local_percentile_500)
                            else float("nan")
                        ),
                        "local_context_windows_500": (
                            score.local_context_windows_500
                        ),
                        "me_support": score.site_support,
                        "overlap_windows": score.overlap_windows,
                        "gc_fraction": gc_fraction,
                        "target_base_fraction": target_fraction,
                        "callable_treatment_count": callable_count,
                        "sampling_weight": record.sampling_weight,
                        "is_case": case is not None,
                        "significant_in_treatment": (
                            case is not None
                            and treatment in case.treatments
                        ),
                    }
                )

    return pd.DataFrame(rows)


def run_continuous_response_models(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if frame.empty:
        return pd.DataFrame()

    for mod_type, mod_frame in frame.groupby("mod_type"):
        for predictor, endpoint in (
            ("site_me", "site_overlapping_ME"),
            ("domain_me", "neighbourhood_500_ME"),
        ):
            result = weighted_cluster_regression(
                mod_frame,
                predictor=predictor,
                include_treatment_fixed_effects=True,
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "endpoint": endpoint,
                    "analysis_scope": "pooled_treatment_adjusted",
                    "treatment": "ALL",
                }
            )
            rows.append(result)

            for treatment, treatment_frame in mod_frame.groupby(
                "treatment"
            ):
                treatment_result = weighted_cluster_regression(
                    treatment_frame,
                    predictor=predictor,
                    include_treatment_fixed_effects=False,
                )
                treatment_result.update(
                    {
                        "mod_type": mod_type,
                        "endpoint": endpoint,
                        "analysis_scope": "treatment_specific",
                        "treatment": treatment,
                    }
                )
                rows.append(treatment_result)

    result_frame = pd.DataFrame(rows)
    if result_frame.empty:
        return result_frame

    for required_column in (
        "cluster_robust_p_two_sided",
        "coefficient_per_0.1_ME",
        "status",
        "analysis_scope",
    ):
        if required_column not in result_frame.columns:
            result_frame[required_column] = np.nan

    result_frame["holm_p_pooled"] = np.nan
    pooled_mask = (
        result_frame["analysis_scope"]
        == "pooled_treatment_adjusted"
    )
    pooled_indices = result_frame.index[pooled_mask].tolist()
    pooled_p = [
        (
            float(result_frame.loc[index, "cluster_robust_p_two_sided"])
            if (
                str(result_frame.loc[index, "status"])
                == "valid_cluster_robust_continuous_model"
                and math.isfinite(float(
                    result_frame.loc[
                        index,
                        "cluster_robust_p_two_sided",
                    ]
                ))
            )
            else 1.0
        )
        for index in pooled_indices
    ]
    for index, adjusted in zip(
        pooled_indices,
        holm_adjust(pooled_p),
    ):
        result_frame.loc[index, "holm_p_pooled"] = adjusted

    result_frame["bh_p_treatment_specific"] = np.nan
    treatment_mask = (
        result_frame["analysis_scope"]
        == "treatment_specific"
    )
    treatment_indices = result_frame.index[
        treatment_mask
    ].tolist()
    treatment_p = [
        (
            float(result_frame.loc[index, "cluster_robust_p_two_sided"])
            if (
                str(result_frame.loc[index, "status"])
                == "valid_cluster_robust_continuous_model"
                and math.isfinite(float(
                    result_frame.loc[
                        index,
                        "cluster_robust_p_two_sided",
                    ]
                ))
            )
            else 1.0
        )
        for index in treatment_indices
    ]
    for index, adjusted in zip(
        treatment_indices,
        bh_adjust(treatment_p),
    ):
        result_frame.loc[
            index,
            "bh_p_treatment_specific",
        ] = adjusted

    interpretations = []
    for _, row in result_frame.iterrows():
        if (
            row.get("status")
            != "valid_cluster_robust_continuous_model"
        ):
            interpretations.append(str(row.get("status")))
            continue
        coefficient = safe_float(
            row.get("coefficient_per_0.1_ME")
        )
        adjusted = (
            safe_float(row.get("holm_p_pooled"))
            if row.get("analysis_scope")
            == "pooled_treatment_adjusted"
            else safe_float(
                row.get("bh_p_treatment_specific")
            )
        )
        if coefficient is None or adjusted is None:
            interpretations.append("not_interpretable")
        elif adjusted < 0.05 and coefficient > 0:
            interpretations.append(
                "higher_baseline_ME_predicts_larger_abs_deltaPL"
            )
        elif adjusted < 0.05 and coefficient < 0:
            interpretations.append(
                "higher_baseline_ME_predicts_smaller_abs_deltaPL"
            )
        else:
            interpretations.append(
                "no_adjusted_continuous_association"
            )
    result_frame["interpretation"] = interpretations
    return result_frame



# =============================================================================
# COMPREHENSIVE INCLUSIVE ASSOCIATION MODELS (VERSION 5)
# =============================================================================

def _normal_or_t_one_sided(
    statistic: float,
    degrees_of_freedom: int,
    direction: str,
) -> float:
    """One-sided probability matching the requested sign."""
    two_sided = two_sided_t_pvalue(statistic, degrees_of_freedom)
    if not math.isfinite(two_sided):
        return float("nan")
    if direction == "greater":
        return two_sided / 2.0 if statistic >= 0 else 1.0 - two_sided / 2.0
    if direction == "less":
        return two_sided / 2.0 if statistic <= 0 else 1.0 - two_sided / 2.0
    raise ValueError(direction)


def _add_fixed_effects(
    columns: list[np.ndarray],
    names: list[str],
    subset: pd.DataFrame,
    column: str,
) -> None:
    levels = sorted(str(value) for value in subset[column].dropna().unique())
    for level in levels[1:]:
        columns.append(
            (subset[column].astype(str).to_numpy() == level).astype(float)
        )
        names.append(f"{column}[{level}]")


def _build_association_design(
    subset: pd.DataFrame,
    predictor: str,
    weights: np.ndarray,
    *,
    adjusted: bool,
    include_treatment_fixed_effects: bool,
) -> tuple[np.ndarray, list[str]]:
    predictor_values = subset[predictor].to_numpy(dtype=float)
    columns = [
        np.ones(len(subset), dtype=float),
        predictor_values,
    ]
    names = ["intercept", predictor]

    if adjusted:
        pl = subset["control_pl"].to_numpy(dtype=float) - 0.5
        log_depth = np.log1p(
            subset["control_depth"].to_numpy(dtype=float)
        )
        log_support = np.log1p(
            subset["me_support"].to_numpy(dtype=float)
        )
        gc = subset["gc_fraction"].to_numpy(dtype=float)
        target = subset["target_base_fraction"].to_numpy(dtype=float)
        overlap = np.log1p(
            subset["overlap_windows"].to_numpy(dtype=float)
        )
        callable_count = subset[
            "callable_treatment_count"
        ].to_numpy(dtype=float)

        nuisance = [
            ("baseline_PL_linear", pl),
            ("baseline_PL_quadratic", pl ** 2),
            ("baseline_PL_cubic", pl ** 3),
            ("log_control_depth", log_depth),
            ("log_ME_support", log_support),
            ("GC_fraction", gc),
            ("target_base_fraction", target),
            ("log_overlap_windows", overlap),
            ("callable_treatment_count", callable_count),
        ]
        for name, values in nuisance:
            columns.append(weighted_center_scale(values, weights))
            names.append(name)

        if include_treatment_fixed_effects and "treatment" in subset.columns:
            _add_fixed_effects(
                columns,
                names,
                subset,
                "treatment",
            )

    return np.column_stack(columns), names


def _association_subset(
    frame: pd.DataFrame,
    *,
    outcome: str,
    predictor: str,
    adjusted: bool,
) -> pd.DataFrame:
    required = [
        outcome,
        predictor,
        "sampling_weight",
        "genomic_block",
    ]
    if adjusted:
        required.extend(
            [
                "control_pl",
                "control_depth",
                "me_support",
                "gc_fraction",
                "target_base_fraction",
                "overlap_windows",
                "callable_treatment_count",
            ]
        )
    subset = frame.dropna(subset=required).copy()
    numeric = [column for column in required if column != "genomic_block"]
    for column in numeric:
        values = pd.to_numeric(subset[column], errors="coerce")
        subset = subset[np.isfinite(values)]
    subset = subset[
        pd.to_numeric(
            subset["sampling_weight"],
            errors="coerce",
        ) > 0
    ]
    return subset


def weighted_cluster_linear_association(
    frame: pd.DataFrame,
    *,
    outcome: str,
    predictor: str,
    adjusted: bool,
    include_treatment_fixed_effects: bool,
) -> dict[str, object]:
    """
    Weighted linear association with genomic-block CR1 covariance.

    This is used for the deliberately inclusive analyses. No case is discarded
    merely because a close matched control does not exist. The adjusted version
    conditions on baseline PL and measured technical/sequence covariates; the
    unadjusted version is a transparent weighted descriptive contrast.
    """
    subset = _association_subset(
        frame,
        outcome=outcome,
        predictor=predictor,
        adjusted=adjusted,
    )
    if len(subset) < 30:
        return {
            "status": "underpowered_fewer_than_30_rows",
            "rows": len(subset),
        }

    y = subset[outcome].to_numpy(dtype=float)
    weights = subset["sampling_weight"].to_numpy(dtype=float)
    weights = weights / float(np.mean(weights))
    X, names = _build_association_design(
        subset,
        predictor,
        weights,
        adjusted=adjusted,
        include_treatment_fixed_effects=include_treatment_fixed_effects,
    )
    n, p = X.shape
    bread = np.linalg.pinv(
        X.T @ (X * weights[:, None]),
        rcond=1e-10,
    )
    beta = bread @ (X.T @ (weights * y))
    residual = y - X @ beta

    labels = subset["genomic_block"].astype(str).to_numpy()
    clusters = np.unique(labels)
    meat = np.zeros((p, p), dtype=float)
    for cluster in clusters:
        mask = labels == cluster
        score = X[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if len(clusters) > 1 and n > p:
        covariance *= (
            len(clusters) / (len(clusters) - 1.0)
        ) * (
            (n - 1.0) / (n - p)
        )

    coefficient = float(beta[1])
    variance = float(covariance[1, 1])
    se = (
        math.sqrt(variance)
        if math.isfinite(variance) and variance >= 0
        else float("nan")
    )
    statistic = coefficient / se if se > 0 else float("nan")
    df = max(len(clusters) - 1, 1)
    p_two = two_sided_t_pvalue(statistic, df)
    p_high = _normal_or_t_one_sided(statistic, df, "greater")
    p_low = _normal_or_t_one_sided(statistic, df, "less")
    critical = t_critical_975(df)
    status = (
        "valid_cluster_robust_linear_model"
        if (
            n >= 100
            and len(clusters) >= 20
            and math.isfinite(se)
            and se > 0
        )
        else "underpowered_or_singular_linear_model"
    )

    return {
        "status": status,
        "model_adjustment": (
            "adjusted_baseline_and_technical"
            if adjusted
            else "unadjusted_weighted_descriptive"
        ),
        "rows": n,
        "weighted_effective_rows": effective_sample_size(weights),
        "genomic_blocks": len(clusters),
        "outcome": outcome,
        "predictor": predictor,
        "coefficient": coefficient,
        "standard_error": se,
        "ci_low": (
            coefficient - critical * se
            if math.isfinite(se)
            else float("nan")
        ),
        "ci_high": (
            coefficient + critical * se
            if math.isfinite(se)
            else float("nan")
        ),
        "cluster_robust_t": statistic,
        "cluster_robust_df": df,
        "cluster_robust_p_high": p_high,
        "cluster_robust_p_low": p_low,
        "cluster_robust_p_two_sided": p_two,
        "outcome_mean": weighted_mean(y, weights),
        "predictor_mean": weighted_mean(
            subset[predictor].to_numpy(dtype=float),
            weights,
        ),
        "model_terms": ",".join(names),
        "inference_scope": (
            "inclusive_callable_genomic_association_not_biological_"
            "replicate_inference"
        ),
    }


def weighted_cluster_logistic_association(
    frame: pd.DataFrame,
    *,
    outcome: str,
    predictor: str,
    adjusted: bool,
    include_treatment_fixed_effects: bool,
    predictor_scale: float,
) -> dict[str, object]:
    """
    Weighted logistic regression fitted by IRLS with genomic-block CR1
    covariance. It estimates whether baseline ME or peak status predicts a
    later threshold-positive methylation response.
    """
    subset = _association_subset(
        frame,
        outcome=outcome,
        predictor=predictor,
        adjusted=adjusted,
    )
    if len(subset) < 30:
        return {
            "status": "underpowered_fewer_than_30_rows",
            "rows": len(subset),
        }

    y = subset[outcome].astype(float).to_numpy()
    event_count = int(np.sum(y > 0.5))
    non_event_count = int(np.sum(y <= 0.5))
    if event_count < 5 or non_event_count < 20:
        return {
            "status": "underpowered_too_few_events_or_nonevents",
            "rows": len(subset),
            "events": event_count,
            "non_events": non_event_count,
        }

    weights = subset["sampling_weight"].to_numpy(dtype=float)
    weights = weights / float(np.mean(weights))
    X, names = _build_association_design(
        subset,
        predictor,
        weights,
        adjusted=adjusted,
        include_treatment_fixed_effects=include_treatment_fixed_effects,
    )
    n, p = X.shape
    beta = np.zeros(p, dtype=float)
    converged = False
    ridge = np.eye(p, dtype=float) * 1e-9
    ridge[0, 0] = 0.0

    for _ in range(100):
        eta = np.clip(X @ beta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        variance = np.maximum(mu * (1.0 - mu), 1e-8)
        information = X.T @ (
            X * (weights * variance)[:, None]
        ) + ridge
        score = X.T @ (weights * (y - mu))
        step = np.linalg.pinv(information, rcond=1e-10) @ score
        beta_new = beta + step
        if float(np.max(np.abs(step))) < 1e-8:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    eta = np.clip(X @ beta, -30.0, 30.0)
    mu = 1.0 / (1.0 + np.exp(-eta))
    variance = np.maximum(mu * (1.0 - mu), 1e-8)
    bread = np.linalg.pinv(
        X.T @ (X * (weights * variance)[:, None]) + ridge,
        rcond=1e-10,
    )
    labels = subset["genomic_block"].astype(str).to_numpy()
    clusters = np.unique(labels)
    meat = np.zeros((p, p), dtype=float)
    for cluster in clusters:
        mask = labels == cluster
        cluster_score = X[mask].T @ (
            weights[mask] * (y[mask] - mu[mask])
        )
        meat += np.outer(cluster_score, cluster_score)
    covariance = bread @ meat @ bread
    if len(clusters) > 1 and n > p:
        covariance *= (
            len(clusters) / (len(clusters) - 1.0)
        ) * (
            (n - 1.0) / (n - p)
        )

    coefficient = float(beta[1])
    variance_beta = float(covariance[1, 1])
    se = (
        math.sqrt(variance_beta)
        if math.isfinite(variance_beta) and variance_beta >= 0
        else float("nan")
    )
    statistic = coefficient / se if se > 0 else float("nan")
    df = max(len(clusters) - 1, 1)
    p_two = two_sided_t_pvalue(statistic, df)
    p_high = _normal_or_t_one_sided(statistic, df, "greater")
    p_low = _normal_or_t_one_sided(statistic, df, "less")
    critical = t_critical_975(df)

    scaled_beta = coefficient * predictor_scale
    scaled_low = (coefficient - critical * se) * predictor_scale
    scaled_high = (coefficient + critical * se) * predictor_scale

    def safe_exp(value: float) -> float:
        return float(math.exp(max(-700.0, min(700.0, value))))

    status = (
        "valid_cluster_robust_logistic_model"
        if (
            converged
            and event_count >= 10
            and non_event_count >= 100
            and len(clusters) >= 20
            and math.isfinite(se)
            and se > 0
        )
        else (
            "exploratory_small_event_logistic_model"
            if (
                converged
                and event_count >= 5
                and len(clusters) >= 10
                and math.isfinite(se)
                and se > 0
            )
            else "underpowered_or_singular_logistic_model"
        )
    )

    return {
        "status": status,
        "model_adjustment": (
            "adjusted_baseline_and_technical"
            if adjusted
            else "unadjusted_weighted_descriptive"
        ),
        "rows": n,
        "events": event_count,
        "non_events": non_event_count,
        "weighted_event_prevalence": weighted_mean(y, weights),
        "weighted_effective_rows": effective_sample_size(weights),
        "genomic_blocks": len(clusters),
        "outcome": outcome,
        "predictor": predictor,
        "predictor_scale": predictor_scale,
        "log_odds_coefficient_per_scale": scaled_beta,
        "odds_ratio_per_scale": safe_exp(scaled_beta),
        "odds_ratio_ci_low": safe_exp(scaled_low),
        "odds_ratio_ci_high": safe_exp(scaled_high),
        "cluster_robust_t": statistic,
        "cluster_robust_df": df,
        "cluster_robust_p_high": p_high,
        "cluster_robust_p_low": p_low,
        "cluster_robust_p_two_sided": p_two,
        "model_terms": ",".join(names),
        "inference_scope": (
            "inclusive_callable_genomic_association_not_biological_"
            "replicate_inference"
        ),
    }


def add_genome_reference_columns(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add genome-background percentiles and means without using case outcomes to
    define the reference. References are calculated from unique noncase sites
    in the sampled callable risk set, with inverse-reservoir weights.

    Priority is modification+contig+strand; sparse strata fall back to
    modification+contig and then modification alone.
    """
    if frame.empty:
        return frame.copy()
    result = frame.copy()
    site_columns = ["mod_type", "chrom", "start", "end", "strand"]
    unique = result.sort_values(site_columns).drop_duplicates(site_columns)
    unique = unique.copy()
    unique["_site_id"] = (
        unique["mod_type"].astype(str)
        + "|"
        + unique["chrom"].astype(str)
        + ":"
        + unique["start"].astype(str)
        + "-"
        + unique["end"].astype(str)
        + "("
        + unique["strand"].astype(str)
        + ")"
    )

    reference_columns = {
        "site_me": (
            "site_genome_percentile",
            "site_genome_background_mean",
            "site_global_peak90",
            "site_global_peak95",
        ),
        "domain_me": (
            "domain_genome_percentile",
            "domain_genome_background_mean",
            "domain_global_peak90",
            "domain_global_peak95",
        ),
    }
    for target in reference_columns:
        for output in reference_columns[target]:
            unique[output] = np.nan

    grouping_levels = [
        [
            "mod_type",
            "chrom",
            "strand",
            "callable_treatment_count",
        ],
        ["mod_type", "chrom", "strand"],
        ["mod_type", "chrom"],
        ["mod_type"],
    ]

    for target, outputs in reference_columns.items():
        percentile_col, mean_col, peak90_col, peak95_col = outputs
        unassigned = pd.Series(True, index=unique.index)

        for grouping in grouping_levels:
            grouped = unique.loc[unassigned].groupby(
                grouping[0] if len(grouping) == 1 else grouping,
                sort=False,
            )
            for _, query_group in grouped:
                first = query_group.iloc[0]
                mask = pd.Series(True, index=unique.index)
                for column in grouping:
                    mask &= unique[column].astype(str) == str(first[column])
                background = unique[
                    mask
                    & (~unique["is_case"].astype(bool))
                    & unique[target].notna()
                ]
                if len(background) < 50:
                    continue
                values = background[target].to_numpy(dtype=float)
                weights = background["sampling_weight"].to_numpy(dtype=float)
                finite = (
                    np.isfinite(values)
                    & np.isfinite(weights)
                    & (weights > 0)
                )
                values = values[finite]
                weights = weights[finite]
                if len(values) < 50:
                    continue
                order = np.argsort(values)
                sorted_values = values[order]
                sorted_weights = weights[order]
                cumulative = np.cumsum(sorted_weights)
                total = float(cumulative[-1])
                background_mean = weighted_mean(values, weights)

                query_indices = query_group.index[unassigned.loc[
                    query_group.index
                ]]
                query_values = unique.loc[
                    query_indices,
                    target,
                ].to_numpy(dtype=float)
                positions = np.searchsorted(
                    sorted_values,
                    query_values,
                    side="right",
                ) - 1
                percentiles = np.zeros(len(query_values), dtype=float)
                valid_positions = positions >= 0
                percentiles[valid_positions] = (
                    cumulative[positions[valid_positions]] / total
                )
                percentiles[~np.isfinite(query_values)] = np.nan

                unique.loc[query_indices, percentile_col] = percentiles
                unique.loc[query_indices, mean_col] = background_mean
                unique.loc[query_indices, peak90_col] = (
                    percentiles >= 0.90
                ).astype(float)
                unique.loc[query_indices, peak95_col] = (
                    percentiles >= 0.95
                ).astype(float)
                unassigned.loc[query_indices] = False

    unique["site_minus_genome_background_mean"] = (
        unique["site_me"] - unique["site_genome_background_mean"]
    )
    unique["domain_minus_genome_background_mean"] = (
        unique["domain_me"] - unique["domain_genome_background_mean"]
    )

    reference_keep = [
        *site_columns,
        "site_genome_percentile",
        "site_genome_background_mean",
        "site_global_peak90",
        "site_global_peak95",
        "domain_genome_percentile",
        "domain_genome_background_mean",
        "domain_global_peak90",
        "domain_global_peak95",
        "site_minus_genome_background_mean",
        "domain_minus_genome_background_mean",
    ]
    result = result.merge(
        unique[reference_keep],
        on=site_columns,
        how="left",
        validate="many_to_one",
    )
    return result


def build_union_site_frame(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse treatment-site rows to one ever-responsive row per site."""
    if frame.empty:
        return frame.copy()
    group_columns = ["mod_type", "chrom", "start", "end", "strand"]
    invariant_first = [
        "context",
        "genomic_block",
        "site_me",
        "domain_me",
        "context_500_me",
        "local_prominence_500",
        "local_percentile_500",
        "local_peak90",
        "local_peak95",
        "local_context_windows_500",
        "me_support",
        "overlap_windows",
        "gc_fraction",
        "target_base_fraction",
        "callable_treatment_count",
        "sampling_weight",
        "site_genome_percentile",
        "site_genome_background_mean",
        "site_global_peak90",
        "site_global_peak95",
        "domain_genome_percentile",
        "domain_genome_background_mean",
        "domain_global_peak90",
        "domain_global_peak95",
        "site_minus_genome_background_mean",
        "domain_minus_genome_background_mean",
    ]
    rows: list[dict[str, object]] = []
    for key, group in frame.groupby(group_columns, sort=False):
        row = dict(zip(group_columns, key))
        for column in invariant_first:
            row[column] = group[column].iloc[0] if column in group else np.nan
        row["control_pl"] = float(group["control_pl"].mean())
        row["control_depth"] = float(group["control_depth"].mean())
        row["treatment_pl"] = float(group["treatment_pl"].mean())
        row["delta_pl"] = float(group["delta_pl"].mean())
        row["abs_delta_pl"] = float(group["abs_delta_pl"].max())
        row["is_case"] = bool(group["is_case"].astype(bool).any())
        row["significant_in_treatment"] = bool(
            group["significant_in_treatment"].astype(bool).any()
        )
        row["significant_treatment_count"] = int(
            group["significant_in_treatment"].astype(bool).sum()
        )
        row["treatment"] = "ANY_TESTED_TREATMENT"
        rows.append(row)
    return pd.DataFrame(rows)


def run_unadjusted_continuous_models(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if frame.empty:
        return pd.DataFrame()
    for mod_type, mod_frame in frame.groupby("mod_type"):
        for predictor, endpoint in (
            ("site_me", "site_overlapping_ME"),
            ("domain_me", "neighbourhood_500_ME"),
        ):
            result = weighted_cluster_linear_association(
                mod_frame,
                outcome="abs_delta_pl",
                predictor=predictor,
                adjusted=False,
                include_treatment_fixed_effects=False,
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "endpoint": endpoint,
                    "analysis_scope": "pooled_treatment_unadjusted",
                    "treatment": "ALL",
                    "coefficient_per_0.1_ME": (
                        float(result.get("coefficient", np.nan)) * 0.1
                    ),
                    "percentage_point_abs_deltaPL_per_0.1_ME": (
                        float(result.get("coefficient", np.nan)) * 10.0
                    ),
                    "cluster_robust_p_two_sided": result.get(
                        "cluster_robust_p_two_sided",
                        np.nan,
                    ),
                }
            )
            rows.append(result)
    output = pd.DataFrame(rows)
    if not output.empty:
        output["holm_p_pooled"] = holm_adjust(
            output["cluster_robust_p_two_sided"].fillna(1.0).tolist()
        )
        output["interpretation"] = [
            (
                "higher_baseline_ME_associated_with_larger_abs_deltaPL"
                if (
                    safe_float(row.get("coefficient")) is not None
                    and float(row.get("coefficient")) > 0
                    and float(row.get("holm_p_pooled", 1.0)) < 0.05
                )
                else (
                    "higher_baseline_ME_associated_with_smaller_abs_deltaPL"
                    if (
                        safe_float(row.get("coefficient")) is not None
                        and float(row.get("coefficient")) < 0
                        and float(row.get("holm_p_pooled", 1.0)) < 0.05
                    )
                    else "no_unadjusted_continuous_association"
                )
            )
            for _, row in output.iterrows()
        ]
    return output


def run_binary_significance_models(
    treatment_frame: pd.DataFrame,
    union_frame: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    predictors = [
        ("site_me", "site_overlapping_ME_per_0.1", 0.1),
        ("domain_me", "neighbourhood_500_ME_per_0.1", 0.1),
        (
            "local_prominence_500",
            "site_minus_own_500bp_context_per_0.1",
            0.1,
        ),
        ("local_peak90", "local_top_decile_peak_yes_vs_no", 1.0),
        ("site_global_peak90", "genome_top_decile_yes_vs_no", 1.0),
    ]
    # The primary threshold-positive estimand is the union of unique sites
    # significant in any tested treatment. Treatment-specific behaviour is
    # already retained in the continuous and historical sensitivity modules.
    for analysis_unit, frame in (
        ("union_site", union_frame),
    ):
        if frame.empty:
            continue
        outcome = (
            "significant_in_treatment"
            if analysis_unit == "treatment_site"
            else "is_case"
        )
        for mod_type, mod_frame in frame.groupby("mod_type"):
            for predictor, label, scale in predictors:
                if predictor not in mod_frame:
                    continue
                for adjusted in (False, True):
                    result = weighted_cluster_logistic_association(
                        mod_frame,
                        outcome=outcome,
                        predictor=predictor,
                        adjusted=adjusted,
                        include_treatment_fixed_effects=(
                            analysis_unit == "treatment_site"
                        ),
                        predictor_scale=scale,
                    )
                    result.update(
                        {
                            "mod_type": mod_type,
                            "analysis_unit": analysis_unit,
                            "hypothesis_predictor": label,
                        }
                    )
                    rows.append(result)

    output = pd.DataFrame(rows)
    if output.empty:
        return output
    output["holm_p_within_design"] = np.nan
    for (_, adjustment), group in output.groupby(
        ["analysis_unit", "model_adjustment"],
    ):
        indices = group.index.tolist()
        pvalues = [
            (
                float(output.loc[index, "cluster_robust_p_two_sided"])
                if (
                    str(output.loc[index, "status"]).startswith("valid_")
                    or str(output.loc[index, "status"]).startswith(
                        "exploratory_"
                    )
                )
                and math.isfinite(float(
                    output.loc[index, "cluster_robust_p_two_sided"]
                ))
                else 1.0
            )
            for index in indices
        ]
        for index, adjusted_p in zip(indices, holm_adjust(pvalues)):
            output.loc[index, "holm_p_within_design"] = adjusted_p

    output["interpretation"] = [
        (
            "higher_predictor_increases_odds_of_later_significant_site"
            if (
                safe_float(row.get("odds_ratio_per_scale")) is not None
                and float(row.get("odds_ratio_per_scale")) > 1
                and float(row.get("holm_p_within_design", 1.0)) < 0.05
            )
            else (
                "higher_predictor_decreases_odds_of_later_significant_site"
                if (
                    safe_float(row.get("odds_ratio_per_scale")) is not None
                    and float(row.get("odds_ratio_per_scale")) < 1
                    and float(row.get("holm_p_within_design", 1.0)) < 0.05
                )
                else "no_multiplicity_adjusted_binary_association"
            )
        )
        for _, row in output.iterrows()
    ]
    return output


def run_case_me_difference_models(
    treatment_frame: pd.DataFrame,
    union_frame: pd.DataFrame,
) -> pd.DataFrame:
    """
    Reverse-orientation companion analysis: quantify how much higher or lower
    baseline ME is at sites that later cross the response threshold.
    """
    rows: list[dict[str, object]] = []
    outcomes = [
        ("site_me", "site_overlapping_ME"),
        ("domain_me", "neighbourhood_500_ME"),
        (
            "local_prominence_500",
            "site_minus_own_500bp_context",
        ),
    ]
    for analysis_unit, frame in (
        ("union_site", union_frame),
    ):
        if frame.empty:
            continue
        predictor = (
            "significant_in_treatment"
            if analysis_unit == "treatment_site"
            else "is_case"
        )
        for mod_type, mod_frame in frame.groupby("mod_type"):
            for outcome, endpoint in outcomes:
                for adjusted in (False, True):
                    result = weighted_cluster_linear_association(
                        mod_frame,
                        outcome=outcome,
                        predictor=predictor,
                        adjusted=adjusted,
                        include_treatment_fixed_effects=(
                            analysis_unit == "treatment_site"
                        ),
                    )
                    result.update(
                        {
                            "mod_type": mod_type,
                            "analysis_unit": analysis_unit,
                            "ME_endpoint": endpoint,
                            "coefficient_interpretation": (
                                "mean_ME_difference_threshold_positive_"
                                "minus_background"
                            ),
                        }
                    )
                    rows.append(result)
    output = pd.DataFrame(rows)
    if output.empty:
        return output
    output["holm_p_within_design"] = np.nan
    for (_, adjustment), group in output.groupby(
        ["analysis_unit", "model_adjustment"],
    ):
        indices = group.index.tolist()
        pvalues = [
            (
                float(output.loc[index, "cluster_robust_p_two_sided"])
                if str(output.loc[index, "status"]).startswith("valid_")
                and math.isfinite(float(
                    output.loc[index, "cluster_robust_p_two_sided"]
                ))
                else 1.0
            )
            for index in indices
        ]
        for index, adjusted_p in zip(indices, holm_adjust(pvalues)):
            output.loc[index, "holm_p_within_design"] = adjusted_p
    output["interpretation"] = [
        (
            "later_significant_sites_have_higher_ME"
            if (
                safe_float(row.get("coefficient")) is not None
                and float(row.get("coefficient")) > 0
                and float(row.get("holm_p_within_design", 1.0)) < 0.05
            )
            else (
                "later_significant_sites_have_lower_ME"
                if (
                    safe_float(row.get("coefficient")) is not None
                    and float(row.get("coefficient")) < 0
                    and float(row.get("holm_p_within_design", 1.0)) < 0.05
                )
                else "no_multiplicity_adjusted_ME_difference"
            )
        )
        for _, row in output.iterrows()
    ]
    return output


def cluster_robust_mean_test(
    frame: pd.DataFrame,
    *,
    value_column: str,
) -> dict[str, object]:
    subset = frame.dropna(
        subset=[value_column, "genomic_block"]
    ).copy()
    values = subset[value_column].to_numpy(dtype=float)
    finite = np.isfinite(values)
    subset = subset.loc[finite].copy()
    values = values[finite]
    if len(values) < 5:
        return {
            "status": "underpowered_fewer_than_5_sites",
            "sites": len(values),
        }

    labels = subset["genomic_block"].astype(str).to_numpy()
    clusters = np.unique(labels)
    mean_value = float(np.mean(values))
    cluster_sums = []
    for cluster in clusters:
        mask = labels == cluster
        cluster_sums.append(float(np.sum(values[mask] - mean_value)))
    # Intercept-only CR1 sandwich variance.
    bread = 1.0 / float(len(values))
    meat = float(np.sum(np.square(cluster_sums)))
    variance = bread * bread * meat
    if len(clusters) > 1:
        variance *= len(clusters) / (len(clusters) - 1.0)
    se = math.sqrt(max(variance, 0.0))
    statistic = mean_value / se if se > 0 else float("nan")
    df = max(len(clusters) - 1, 1)
    critical = t_critical_975(df)
    status = (
        "valid_cluster_robust_direct_contrast"
        if len(values) >= 10 and len(clusters) >= 5 and se > 0
        else "exploratory_small_direct_contrast"
    )
    return {
        "status": status,
        "sites": len(values),
        "genomic_blocks": len(clusters),
        "value_column": value_column,
        "mean_contrast": mean_value,
        "median_contrast": float(np.median(values)),
        "fraction_positive": float(np.mean(values > 0)),
        "standard_error": se,
        "ci_low": mean_value - critical * se,
        "ci_high": mean_value + critical * se,
        "cluster_robust_t": statistic,
        "cluster_robust_df": df,
        "p_high": _normal_or_t_one_sided(statistic, df, "greater"),
        "p_low": _normal_or_t_one_sided(statistic, df, "less"),
        "p_two_sided": two_sided_t_pvalue(statistic, df),
        "inference_scope": (
            "within_candidate_genomic_location_contrast_not_biological_"
            "replicate_inference"
        ),
    }


def run_direct_peak_contrasts(
    union_frame: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if union_frame.empty:
        return pd.DataFrame()
    cases = union_frame[union_frame["is_case"].astype(bool)].copy()
    for mod_type, mod_frame in cases.groupby("mod_type"):
        for column, hypothesis, plain in (
            (
                "local_prominence_500",
                "site_overlapping_ME_vs_own_500bp_context",
                (
                    "Is ME in windows physically overlapping the later "
                    "responsive base higher than non-overlapping windows "
                    "within +/-500 bp?"
                ),
            ),
            (
                "site_minus_genome_background_mean",
                "site_overlapping_ME_vs_callable_genome_background",
                (
                    "Is ME at the later responsive base higher than the "
                    "eligible callable genome-wide background?"
                ),
            ),
        ):
            result = cluster_robust_mean_test(
                mod_frame,
                value_column=column,
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "hypothesis": hypothesis,
                    "plain_language_question": plain,
                    "median_local_percentile_500": (
                        float(mod_frame["local_percentile_500"].median())
                        if "local_percentile_500" in mod_frame
                        else float("nan")
                    ),
                    "fraction_local_peak90": (
                        float(mod_frame["local_peak90"].mean())
                        if "local_peak90" in mod_frame
                        else float("nan")
                    ),
                    "fraction_local_peak95": (
                        float(mod_frame["local_peak95"].mean())
                        if "local_peak95" in mod_frame
                        else float("nan")
                    ),
                    "median_genome_percentile": (
                        float(mod_frame["site_genome_percentile"].median())
                        if "site_genome_percentile" in mod_frame
                        else float("nan")
                    ),
                    "fraction_genome_peak90": (
                        float(mod_frame["site_global_peak90"].mean())
                        if "site_global_peak90" in mod_frame
                        else float("nan")
                    ),
                    "fraction_genome_peak95": (
                        float(mod_frame["site_global_peak95"].mean())
                        if "site_global_peak95" in mod_frame
                        else float("nan")
                    ),
                }
            )
            rows.append(result)
    output = pd.DataFrame(rows)
    if not output.empty:
        output["holm_p_high"] = holm_adjust(
            output["p_high"].fillna(1.0).tolist()
        )
        output["holm_p_two_sided"] = holm_adjust(
            output["p_two_sided"].fillna(1.0).tolist()
        )
        output["interpretation"] = [
            (
                "candidate_site_ME_higher_than_reference"
                if (
                    safe_float(row.get("mean_contrast")) is not None
                    and float(row.get("mean_contrast")) > 0
                    and float(row.get("holm_p_high", 1.0)) < 0.05
                )
                else (
                    "candidate_site_ME_lower_than_reference"
                    if (
                        safe_float(row.get("mean_contrast")) is not None
                        and float(row.get("mean_contrast")) < 0
                        and float(row.get("holm_p_two_sided", 1.0)) < 0.05
                    )
                    else "no_multiplicity_adjusted_direct_contrast"
                )
            )
            for _, row in output.iterrows()
        ]
    return output


def write_hypothesis_results_report(
    output_dir: Path,
    *,
    strict_matched: pd.DataFrame,
    continuous: pd.DataFrame,
    binary: pd.DataFrame,
    me_differences: pd.DataFrame,
    direct: pd.DataFrame,
    union_sites: pd.DataFrame,
) -> None:
    """
    Human-readable report organized by biological hypothesis rather than by
    implementation detail. It is intentionally suitable for upload and rapid
    interpretation after a run.
    """
    def f(value: object, digits: int = 4) -> str:
        number = safe_float(value)
        if number is None:
            return "NA"
        if abs(number) < 0.0001 and number != 0:
            return f"{number:.3e}"
        return f"{number:.{digits}f}"

    lines = [
        "ME-NODE ENRICHMENT VERSION 5: HYPOTHESES AND OUTCOMES",
        "=" * 64,
        "",
        "Terminology",
        "-----------",
        "Site-overlapping ME: support-weighted median ME from true Modkit",
        "windows that physically overlap the exact base.",
        "Local +/-500-bp context: median ME of non-overlapping windows whose",
        "midpoints lie within 500 bp of that base.",
        "Genome background: eligible callable sites sampled across the same",
        "modification/contig/strand risk set.",
        "",
        "Adjusted analyses control baseline PL, depth, ME support, local",
        "sequence composition, overlapping-window count, callability and,",
        "where relevant, treatment. Unadjusted analyses retain only sampling",
        "weights and genomic-block uncertainty and are descriptive.",
        "",
    ]

    hypotheses = [
        (
            "H1",
            "Does higher pre-treatment site-overlapping ME predict a larger later |delta PL|?",
            "continuous",
            "site_overlapping_ME",
        ),
        (
            "H2",
            "Does higher pre-treatment +/-500-bp neighbourhood ME predict a larger later |delta PL|?",
            "continuous",
            "neighbourhood_500_ME",
        ),
        (
            "H3",
            "Does pre-treatment site-overlapping ME predict whether a site later crosses the response threshold?",
            "binary",
            "site_overlapping_ME_per_0.1",
        ),
        (
            "H4",
            "Does pre-treatment +/-500-bp neighbourhood ME predict whether a site later crosses the response threshold?",
            "binary",
            "neighbourhood_500_ME_per_0.1",
        ),
        (
            "H5",
            "How much higher is site-overlapping ME at later threshold-positive sites than at the callable background?",
            "me_difference",
            "site_overlapping_ME",
        ),
        (
            "H6",
            "How much higher is +/-500-bp neighbourhood ME at later threshold-positive sites than at the callable background?",
            "me_difference",
            "neighbourhood_500_ME",
        ),
        (
            "H7",
            "Are later threshold-positive sites focal ME peaks relative to their own +/-500-bp context?",
            "direct",
            "site_overlapping_ME_vs_own_500bp_context",
        ),
        (
            "H8",
            "Are later threshold-positive sites high-ME positions relative to the callable genome background?",
            "direct",
            "site_overlapping_ME_vs_callable_genome_background",
        ),
        (
            "H9",
            "After strict case-specific PL/depth/support matching, do threshold-positive sites remain ME-enriched?",
            "strict",
            "",
        ),
    ]

    for code, question, source, key in hypotheses:
        lines.extend([f"{code}. {question}", "-" * 64])
        if source == "continuous":
            subset = continuous[
                continuous["endpoint"].astype(str) == key
            ] if not continuous.empty else pd.DataFrame()
            if subset.empty:
                lines.append("No testable result.")
            else:
                for _, row in subset.sort_values(
                    ["mod_type", "analysis_scope"]
                ).iterrows():
                    lines.append(
                        f"{row.get('mod_type')} | "
                        f"{row.get('analysis_scope')}: "
                        f"effect per 0.10 ME={f(row.get('coefficient_per_0.1_ME'))}, "
                        f"p={f(row.get('cluster_robust_p_two_sided'))}, "
                        f"adjusted p={f(row.get('holm_p_pooled'))}, "
                        f"status={row.get('status')}"
                    )
        elif source == "binary":
            subset = binary[
                binary["hypothesis_predictor"].astype(str) == key
            ] if not binary.empty else pd.DataFrame()
            if subset.empty:
                lines.append("No testable result.")
            else:
                for _, row in subset.sort_values(
                    ["analysis_unit", "mod_type", "model_adjustment"]
                ).iterrows():
                    lines.append(
                        f"{row.get('analysis_unit')} | {row.get('mod_type')} | "
                        f"{row.get('model_adjustment')}: "
                        f"OR={f(row.get('odds_ratio_per_scale'))} "
                        f"[{f(row.get('odds_ratio_ci_low'))}, "
                        f"{f(row.get('odds_ratio_ci_high'))}], "
                        f"p={f(row.get('cluster_robust_p_two_sided'))}, "
                        f"Holm p={f(row.get('holm_p_within_design'))}, "
                        f"events={row.get('events')}, status={row.get('status')}"
                    )
        elif source == "me_difference":
            subset = me_differences[
                me_differences["ME_endpoint"].astype(str) == key
            ] if not me_differences.empty else pd.DataFrame()
            if subset.empty:
                lines.append("No testable result.")
            else:
                for _, row in subset.sort_values(
                    ["analysis_unit", "mod_type", "model_adjustment"]
                ).iterrows():
                    lines.append(
                        f"{row.get('analysis_unit')} | {row.get('mod_type')} | "
                        f"{row.get('model_adjustment')}: "
                        f"threshold-positive minus background ME="
                        f"{f(row.get('coefficient'))} "
                        f"[{f(row.get('ci_low'))}, {f(row.get('ci_high'))}], "
                        f"p={f(row.get('cluster_robust_p_two_sided'))}, "
                        f"Holm p={f(row.get('holm_p_within_design'))}, "
                        f"status={row.get('status')}"
                    )
        elif source == "direct":
            subset = direct[
                direct["hypothesis"].astype(str) == key
            ] if not direct.empty else pd.DataFrame()
            if subset.empty:
                lines.append("No testable result.")
            else:
                for _, row in subset.sort_values("mod_type").iterrows():
                    lines.append(
                        f"{row.get('mod_type')}: mean contrast="
                        f"{f(row.get('mean_contrast'))} "
                        f"[{f(row.get('ci_low'))}, {f(row.get('ci_high'))}], "
                        f"Holm one-sided p={f(row.get('holm_p_high'))}, "
                        f"median local percentile="
                        f"{f(row.get('median_local_percentile_500'))}, "
                        f"local top-decile fraction="
                        f"{f(row.get('fraction_local_peak90'))}, "
                        f"genome top-decile fraction="
                        f"{f(row.get('fraction_genome_peak90'))}, "
                        f"status={row.get('status')}"
                    )
        else:
            if strict_matched.empty:
                lines.append("No strict matched result.")
            else:
                for _, row in strict_matched.sort_values(
                    ["mod_type", "endpoint"]
                ).iterrows():
                    lines.append(
                        f"{row.get('mod_type')} | {row.get('analysis')}: "
                        f"cases={row.get('site_cases')}, "
                        f"effect={f(row.get('observed_mean_effect'))}, "
                        f"Holm p={f(row.get('holm_primary_p'))}, "
                        f"status={row.get('inference_validity')}"
                    )
        lines.append("")

    lines.extend(
        [
            "Descriptive case retention and zero structure",
            "-" * 64,
        ]
    )
    if not union_sites.empty:
        for mod_type, group in union_sites.groupby("mod_type"):
            cases = group[group["is_case"].astype(bool)]
            lines.append(
                f"{mod_type}: scored cases={len(cases)}; "
                f"mean site ME={f(cases['site_me'].mean())}; "
                f"fraction site ME=0={f((cases['site_me'] == 0).mean())}; "
                f"mean +/-500-bp ME={f(cases['domain_me'].mean())}; "
                f"fraction +/-500-bp ME=0="
                f"{f((cases['domain_me'] == 0).mean())}."
            )

    lines.extend(
        [
            "",
            "Interpretive hierarchy",
            "-" * 64,
            "1. Give greatest weight to adjusted inclusive models that remain",
            "   significant after Holm correction and show stable direction.",
            "2. Use unadjusted results to show the raw genomic pattern, not as",
            "   proof that ME acts independently of baseline PL/depth.",
            "3. Strict matched results are the strongest case-control design",
            "   when adequately powered; underpowered strict matching does not",
            "   erase an inclusive adjusted association.",
            "4. These are genomic-location statistics from pooled control ME.",
            "   They do not establish biological-replicate reproducibility,",
            "   stable epialleles or coherent whole-genome subpopulations.",
            "",
        ]
    )

    (output_dir / "HYPOTHESES_AND_RESULTS.txt").write_text(
        "\n".join(lines) + "\n"
    )

# =============================================================================
# PLOTS AND REPORTS
# =============================================================================

def plot_case_percentiles(
    case_frame: pd.DataFrame,
    output: Path,
) -> None:
    if case_frame.empty:
        return
    for endpoint, column in (
        ("local", "local_percentile"),
        ("global", "global_percentile"),
    ):
        values = [
            group[column].dropna().to_numpy(dtype=float)
            for _, group in case_frame.groupby("mod_type")
        ]
        labels = [
            str(label)
            for label, _ in case_frame.groupby("mod_type")
        ]
        if not values or not any(len(item) for item in values):
            continue
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        try:
            ax.boxplot(values, tick_labels=labels)
        except TypeError:
            # Matplotlib <3.9 uses ``labels``.
            ax.boxplot(values, labels=labels)
        ax.axhline(0.5, linewidth=0.8, linestyle="--")
        ax.axhline(0.9, linewidth=0.8, linestyle=":")
        ax.set_ylabel(
            "Case percentile within eligible "
            + (
                "local neighbourhood"
                if endpoint == "local"
                else "genome-wide matched background"
            )
        )
        ax.set_xlabel("Modification label")
        fig.tight_layout()
        fig.savefig(
            output.parent / f"{endpoint}_case_percentiles.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)


def interpretation_label(
    effect: float,
    adjusted_high: float,
    adjusted_low: float,
    adjusted_two_sided: float,
) -> str:
    if not all(
        math.isfinite(value)
        for value in (
            effect,
            adjusted_high,
            adjusted_low,
            adjusted_two_sided,
        )
    ):
        return "untestable"
    if effect > 0 and adjusted_high < 0.05:
        return "supports_higher_control_ME"
    if (
        effect < 0
        and adjusted_low < 0.05
        and adjusted_two_sided < 0.05
    ):
        return "significant_opposite_direction_lower_control_ME"
    return "no_detectable_enrichment"


def write_summary(
    output_dir: Path,
    primary: pd.DataFrame,
    case_qc: pd.DataFrame,
    me_qc: pd.DataFrame,
    metadata_qc: pd.DataFrame,
    me_inputs: Sequence[tuple[str, str, Path]],
    provenance_qc: pd.DataFrame,
    direction_results: Optional[pd.DataFrame] = None,
    treatment_results: Optional[pd.DataFrame] = None,
    sensitivity_results: Optional[pd.DataFrame] = None,
    continuous_results: Optional[pd.DataFrame] = None,
) -> None:
    lines = [
        "ME-node enrichment version 5",
        "============================",
        "",
        "Questions",
        "---------",
        "1. Local specificity: are responsive sites higher in control ME than",
        "   other eligible callable sites in the same neighbourhood?",
        "2. Regional enrichment: are responsive sites located in higher-control-ME",
        "   domains than random eligible callable sites elsewhere in the genome?",
        "",
        "Inference scope",
        "---------------",
        "P-values are case-specific matched-set genomic randomization results for",
        "the supplied control-ME tracks. Genomic sites are not biological replicates.",
        "With one pooled control track per modification, this analysis does not",
        "estimate control-to-control biological variability or treatment-induced ΔME.",
        "",
        "Primary results",
        "---------------",
    ]

    if primary.empty:
        lines.append("No primary test had enough eligible cases and controls.")
    else:
        for _, row in primary.iterrows():
            endpoint_name = (
                "local site-specificity"
                if row["endpoint"] == "local"
                else "global domain enrichment"
            )
            if (
                row.get("inference_validity", "")
                != "valid_matched_set_genomic_test"
            ):
                lines.append(
                    f"{row['mod_type']} — {endpoint_name}: "
                    f"{row.get('inference_validity', 'untestable')}; "
                    f"descriptive effect="
                    f"{float(row.get('observed_mean_effect', float('nan'))):.6g}; "
                    f"uninterpreted one-sided p="
                    f"{float(row.get('genomic_randomization_p_high', float('nan'))):.6g}; "
                    f"max |SMD|="
                    f"{float(row.get('max_abs_covariate_smd', float('nan'))):.3g}; "
                    f"{int(row.get('site_cases', 0))} retained from "
                    f"{int(row.get('site_cases_before_balance_trim', row.get('site_cases', 0)))} "
                    f"matched cases / "
                    f"{int(row.get('event_clusters', 0))} event clusters."
                )
                continue
            lines.extend(
                [
                    (
                        f"{row['mod_type']} — {endpoint_name}: "
                        f"effect={row['observed_mean_effect']:.6g}; "
                        f"95% cluster-bootstrap CI "
                        f"{row['bootstrap_ci_low']:.6g} to "
                        f"{row['bootstrap_ci_high']:.6g}; "
                        f"one-sided p={row['genomic_randomization_p_high']:.6g}; "
                        f"Holm p={row['holm_primary_p_high']:.6g}; "
                        f"max |SMD|={row.get('max_abs_covariate_smd', float('nan')):.3g}; "
                        f"{int(row['site_cases'])} sites retained from "
                        f"{int(row.get('site_cases_before_balance_trim', row['site_cases']))} "
                        f"matched cases / "
                        f"{int(row['event_clusters'])} event clusters."
                    ),
                    f"Interpretation: {row['interpretation']}.",
                ]
            )

    lines.extend(
        [
            "",
            "Model projection from the two primary questions",
            "-----------------------------------------------",
        ]
    )
    if primary.empty:
        lines.append(
            "No modification label had enough valid cases for a model-level "
            "projection."
        )
    else:
        for mod_type, subset in primary.groupby("mod_type"):
            rows = {
                str(row["endpoint"]): row
                for _, row in subset.iterrows()
            }
            local = rows.get("local")
            global_row = rows.get("global")
            if local is None or global_row is None:
                lines.append(
                    f"{mod_type}: incomplete primary pair; interpretation "
                    "requires both local and global tests."
                )
                continue
            local_interpretation = str(local["interpretation"])
            global_interpretation = str(global_row["interpretation"])
            if (
                local_interpretation == "supports_higher_control_ME"
                and global_interpretation == "supports_higher_control_ME"
            ):
                statement = (
                    "responsive sites occupy both unusually high-ME local "
                    "positions and unusually high-ME genomic domains"
                )
            elif (
                local_interpretation == "supports_higher_control_ME"
                and global_interpretation
                == "no_detectable_enrichment"
            ):
                statement = (
                    "responsive sites are locally special within otherwise "
                    "ordinary-ME genomic domains"
                )
            elif (
                local_interpretation == "no_detectable_enrichment"
                and global_interpretation
                == "supports_higher_control_ME"
            ):
                statement = (
                    "responsive sites preferentially occur in high-ME domains, "
                    "but the exact sites are not locally exceptional"
                )
            elif (
                local_interpretation == "no_detectable_enrichment"
                and global_interpretation
                == "no_detectable_enrichment"
            ):
                statement = (
                    "control ME does not detectably predict responsive sites "
                    "at either the local-site or broad-domain level"
                )
            elif "opposite_direction" in (
                local_interpretation + global_interpretation
            ):
                statement = (
                    "at least one primary endpoint is significantly lower "
                    "than the eligible background, opposing the high-ME model"
                )
            else:
                statement = (
                    "the primary pair is underpowered or otherwise "
                    "non-interpretable"
                )
            lines.append(f"{mod_type}: {statement}.")

    lines.extend(
        [
            "",
            "Continuous response models",
            "--------------------------",
            "These models use the random callable risk-set sample and test whether",
            "baseline ME predicts treatment-associated |delta PL| after adjusting",
            "for baseline PL (cubic), depth, ME support, sequence composition,",
            "window count, callability and treatment. Reported coefficients are",
            "changes in |delta PL| per 0.10 increase in baseline ME.",
        ]
    )
    if continuous_results is None or continuous_results.empty:
        lines.append("No continuous response model was estimable.")
    else:
        pooled = continuous_results[
            continuous_results["analysis_scope"]
            == "pooled_treatment_adjusted"
        ]
        for _, row in pooled.iterrows():
            status = str(row.get("status", "untestable"))
            if status != "valid_cluster_robust_continuous_model":
                lines.append(
                    f"{row.get('mod_type', '')}/{row.get('endpoint', '')}: "
                    f"{status}."
                )
                continue
            coefficient = float(row["coefficient_per_0.1_ME"])
            low = float(row["ci_low_per_0.1_ME"])
            high = float(row["ci_high_per_0.1_ME"])
            raw_p = float(row["cluster_robust_p_two_sided"])
            adjusted_p = float(row["holm_p_pooled"])
            lines.append(
                f"{row['mod_type']} — {row['endpoint']}: "
                f"adjusted slope={coefficient:.6g} |delta PL per 0.10 ME "
                f"({coefficient * 100.0:.3g} percentage points); "
                f"95% CI {low:.6g} to {high:.6g}; "
                f"cluster-robust p={raw_p:.6g}; Holm p={adjusted_p:.6g}; "
                f"partial R2={float(row['partial_R2_baseline_ME']):.4g}; "
                f"{int(row['rows'])} site-treatment rows / "
                f"{int(row['genomic_blocks'])} genomic blocks. "
                f"Interpretation: {row['interpretation']}."
            )

    # Summarize exploratory findings without allowing them to overwrite the
    # prespecified primary conclusion.
    lines.extend(
        [
            "",
            "Notable exploratory and robustness results",
            "------------------------------------------",
        ]
    )
    notable: list[str] = []

    def collect_notable(
        frame: Optional[pd.DataFrame],
        label: str,
    ) -> None:
        if frame is None or frame.empty:
            return
        for _, row in frame.iterrows():
            if (
                row.get("inference_validity", "")
                not in {"", "valid_matched_set_genomic_test"}
            ):
                continue
            effect = safe_float(row.get("observed_mean_effect"))
            if effect is None:
                continue
            high_adj = safe_float(row.get("bh_high_p"))
            low_adj = safe_float(row.get("bh_low_p"))
            two_adj = safe_float(row.get("bh_two_sided_p"))
            high_raw = safe_float(
                row.get("genomic_randomization_p_high")
            )
            low_raw = safe_float(
                row.get("genomic_randomization_p_low")
            )
            descriptor_parts = [
                str(row.get("mod_type", "")),
                str(row.get("endpoint", "")),
            ]
            for field in (
                "treatment",
                "direction",
                "sensitivity",
                "sensitivity_value",
            ):
                value = row.get(field, "")
                if str(value) not in {"", "nan", "None"}:
                    descriptor_parts.append(f"{field}={value}")
            descriptor = ", ".join(
                part for part in descriptor_parts if part
            )

            if (
                effect > 0
                and high_adj is not None
                and high_adj < 0.05
            ):
                notable.append(
                    f"{label}: {descriptor}: higher control ME "
                    f"(effect={effect:.6g}, BH p={high_adj:.6g})."
                )
            elif (
                effect < 0
                and low_adj is not None
                and two_adj is not None
                and low_adj < 0.05
                and two_adj < 0.05
            ):
                notable.append(
                    f"{label}: {descriptor}: significantly lower control ME, "
                    f"opposite to the proposed enrichment "
                    f"(effect={effect:.6g}, BH two-sided p={two_adj:.6g})."
                )
            elif (
                (effect > 0 and high_raw is not None and high_raw < 0.05)
                or (
                    effect < 0
                    and low_raw is not None
                    and low_raw < 0.05
                )
            ):
                notable.append(
                    f"{label}: {descriptor}: nominal unadjusted signal only "
                    f"(effect={effect:.6g}); treat as suggestive."
                )

    collect_notable(direction_results, "direction")
    collect_notable(treatment_results, "treatment")
    collect_notable(sensitivity_results, "sensitivity")

    if notable:
        lines.extend(notable[:15])
        if len(notable) > 15:
            lines.append(
                f"{len(notable) - 15} additional notable rows are retained "
                "in the result tables."
            )
    else:
        lines.append(
            "No exploratory result survived its reported multiple-testing "
            "correction."
        )

    if sensitivity_results is not None and not sensitivity_results.empty:
        valid = sensitivity_results[
            sensitivity_results["inference_validity"]
            == "valid_matched_set_genomic_test"
        ]
        if not valid.empty:
            lines.append("")
            lines.append("Sensitivity direction summary:")
            for (mod_type, endpoint), subset in valid.groupby(
                ["mod_type", "endpoint"]
            ):
                positive = int(
                    (subset["observed_mean_effect"] > 0).sum()
                )
                negative = int(
                    (subset["observed_mean_effect"] < 0).sum()
                )
                significant_high = int(
                    (subset.get("bh_high_p", pd.Series(dtype=float)) < 0.05)
                    .fillna(False)
                    .sum()
                )
                lines.append(
                    f"- {mod_type}/{endpoint}: {positive} positive and "
                    f"{negative} negative estimates across {len(subset)} "
                    f"valid sensitivities; {significant_high} BH-significant "
                    "in the predicted high-ME direction."
                )

    lines.extend(
        [
            "",
            "Case retention",
            "--------------",
        ]
    )
    if not case_qc.empty:
        grouped = case_qc.groupby("mod_type", as_index=False).agg(
            unique_cases=("unique_case", "sum"),
            callable_cases=("callable", "sum"),
            scored_cases=("scored", "sum"),
            analysed_cases=("analysed", "sum"),
        )
        for _, row in grouped.iterrows():
            lines.append(
                f"{row['mod_type']}: unique={int(row['unique_cases'])}; "
                f"callable={int(row['callable_cases'])}; "
                f"window-scored={int(row['scored_cases'])}; "
                f"analysed={int(row['analysed_cases'])}."
            )

    unstranded = me_qc[
        (me_qc.get("format") == "windows.bedgraph")
        & (me_qc.get("strand") == ".")
    ] if not me_qc.empty else pd.DataFrame()
    lines.extend(
        [
            "",
            "Critical interpretation safeguards",
            "----------------------------------",
            (
                "- Standard five-column windows.bedgraph is unstranded. "
                "Accordingly, the analysis is labelled unstranded and does not "
                "claim same-strand specificity."
                if not unstranded.empty
                else "- Stranded custom window tracks were detected."
            ),
            "- A site score uses only entropy windows physically overlapping the",
            "  candidate coordinate. It does not prove that the candidate base was",
            "  one of the positions used by Modkit to define the entropy pattern.",
            "- regions.bed was used only for QC, never as a peak window.",
            "- Background means callable sites absent from the supplied complete",
            "  metadata case union; absence is not independently verified as a",
            "  DMR q-value because the supplied DMR schema has no q-value field.",
            "- Ordinary ME can be influenced by marginal methylation probabilities.",
            "  Read-pattern counts are required to prove coordinated epialleles or",
            "  discrete population substructures.",
            "- Local target-base density is matched as a technical covariate, but",
            "  methyltransferase motif identity/spacing is not inferred without",
            "  explicit upstream motif provenance or a dedicated motif input.",
            (
                "- Upstream Dorado/Modkit provenance was found; inspect "
                "qc/upstream_provenance.tsv before chemistry-specific claims."
                if (
                    not provenance_qc.empty
                    and (
                        provenance_qc["provenance_json_present"].any()
                        or provenance_qc["modkit_log_present"].any()
                    )
                )
                else (
                    "- No machine-readable Dorado/Modkit provenance was found. "
                    "Modification labels are treated as user-supplied labels, "
                    "not independently verified chemistry-specific measurements."
                )
            ),
            "- If the same control reads contributed both to DMR case discovery",
            "  and to control ME, selection and entropy are statistically coupled.",
            "  Independent validation or replicate cross-fitting is preferable.",
            "- The case set is threshold-selected. The test distinguishes supplied",
            "  responsive calls from callable sites absent from that supplied set;",
            "  it does not estimate a latent error-free biological response class.",
            "- Bootstrap intervals resample retained spatial event clusters while",
            "  holding the constructed matched sets fixed. The matched-set genomic",
            "  randomization p-value is the primary inferential quantity.",
            "",
            "Output guide",
            "------------",
            "results/primary_tests.tsv: local and global confirmatory tests",
            "results/continuous_response_models.tsv: adjusted |delta PL| models",
            "results/direction_tests.tsv: methylation gain/loss analyses",
            "results/treatment_tests.tsv: treatment-specific exploratory tests",
            "results/sensitivity_tests.tsv: support, scale and clustering checks",
            "data/case_level_results.tsv.gz: every initially matched case",
            "data/matched_sets.tsv.gz: exact matched-set membership",
            "data/continuous_response_sites.tsv.gz: continuous-model input rows",
            "qc/primary_posttrim_balance.tsv: final primary balance diagnostics",
            "qc/: input, coverage, coordinate and control-pool diagnostics",
        ]
    )

    (output_dir / "SUMMARY.txt").write_text(
        "\n".join(lines) + "\n"
    )


def write_methods_report(
    output_dir: Path,
    args: argparse.Namespace,
    me_qc: pd.DataFrame,
) -> None:
    """Write a concise, manuscript-facing description of the implemented test."""
    strand_mode = (
        "unstranded"
        if (
            not me_qc.empty
            and (me_qc.get("strand") == ".").any()
        )
        else "stranded"
    )
    lines = [
        "ME-node enrichment version 5: analysis methods",
        "================================================",
        "",
        "Case definition",
        "---------------",
        "Cases were the union of exact methylation-change sites reported in at",
        "least one metadata treatment/modification column. Exact duplicates were",
        "collapsed by modification label, contig, start, end and strand. Treatment",
        "recurrence was retained as metadata but was not required for case status.",
        "Bracketed metadata values were interpreted as [control PL, treatment PL],",
        "and direction was calculated as treatment PL minus control PL.",
        "",
        "Entropy inputs",
        "--------------",
        "All formal endpoints used individual windows from windows.bedgraph.",
        "regions.bed was used only for requested/finite coverage and successful/",
        "failed-window QC. Duplicate window intervals were rejected rather than",
        "averaged. Windows below the specified support threshold were excluded.",
        f"The supplied window tracks were treated as {strand_mode}.",
        "",
        "Primary estimands",
        "------------------",
        f"Local specificity: site-associated control ME was compared with eligible",
        f"callable noncase sites within +/- {args.local_radius} bp on the same",
        "contig. Site-associated ME was the support-weighted median of entropy",
        "windows physically overlapping the base.",
        "",
        f"Global domain enrichment: the support-weighted median control ME within",
        f"+/- {args.domain_radius} bp of each case was compared with randomly",
        "sampled eligible callable noncase sites elsewhere on the same contig.",
        "",
        "Background and adjustment",
        "-------------------------",
        "A background site had to be absent from the complete supplied case union,",
        "have the same modification label, contig, DMR strand and canonical base,",
        "and be callable in every treatment in which the corresponding case was",
        "reported. Nearby noncase sites were retained in the local comparison.",
        "Each case was matched with replacement to a small case-specific set of",
        "nearest controls passing hard calipers on control PL, log control depth,",
        "log ME support, GC fraction and target-base fraction. Controls were chosen",
        "without using ME outcomes. Residual aggregate imbalance was handled by",
        "outcome-blind common-support trimming, capped at 50%; formal interpretation",
        "required maximum absolute post-match SMD <= 0.10.",
        "",
        "Inference",
        "---------",
        "Primary p-values came from exact label randomization within each matched",
        "set: the case label was reassigned among the observed case and its matched",
        "controls. Nearby cases were aggregated into spatial event clusters before",
        "forming the mean statistic. Confidence intervals used a fixed-matched-set",
        "cluster bootstrap. High-ME, low-ME and two-sided p-values were reported;",
        "the high-ME direction was prespecified. Holm correction covered the",
        "primary modification-by-local/global family.",
        "",
        "Continuous response analysis",
        "----------------------------",
        "A secondary weighted model used cases plus the uniformly sampled callable",
        "risk-set reservoir to test whether baseline ME predicted treatment-specific",
        "|delta PL|. Weighted least squares adjusted for a cubic baseline-PL term,",
        "depth, ME support, sequence composition, overlap-window count, callability",
        "and treatment. Genomic-block CR1 standard errors were used. This model",
        "avoids defining the outcome solely by a significance threshold.",
        "",
        "Inclusive threshold-positive and direct peak analyses",
        "----------------------------------------------------",
        "Version 5 additionally retained every threshold-positive site with valid ME",
        "rather than requiring a strict matched control. Weighted logistic models",
        "tested whether site-overlapping ME, +/-500-bp neighbourhood ME, local",
        "prominence, local top-decile peak status and genome top-decile status",
        "predicted later threshold-positive response. Weighted linear companion",
        "models quantified the ME difference between threshold-positive and",
        "background sites. Each analysis was run both without PL/depth adjustment",
        "and with the full baseline/technical adjustment set. Direct candidate-only",
        "contrasts compared site-overlapping ME with non-overlapping windows within",
        "+/-500 bp and with the eligible callable genome background. Holm correction",
        "was applied within each prespecified design family.",
        "",
        "Scale definitions",
        "-----------------",
        "Site-overlapping ME uses only true entropy windows physically overlapping",
        "the exact base. The local context is the non-overlapping window distribution",
        "within +/-500 bp. The genome reference is the eligible callable background",
        "on the same modification/contig/strand stratum. A +/-500-bp neighbourhood",
        "ME endpoint is a local domain summary and is not a whole-genome methylotype.",
        "",
        "Scope",
        "-----",
        "These are conditional genomic-location tests for the supplied control-ME",
        "tracks. They do not treat genomic sites as biological replicates, estimate",
        "control-to-control biological variance from a pooled track, test treatment",
        "minus control ME, or establish coordinated epialleles from entropy alone.",
    ]
    (output_dir / "METHODS.txt").write_text("\\n".join(lines) + "\\n")


def zip_results(output_dir: Path) -> Path:
    zip_path = output_dir.with_suffix(".zip")
    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
    ) as archive:
        for path in sorted(output_dir.rglob("*")):
            if path.is_file():
                archive.write(
                    path,
                    arcname=str(path.relative_to(output_dir.parent)),
                )
    return zip_path


# =============================================================================
# MAIN
# =============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    if args.min_window_support < 1:
        raise ValueError("--min-window-support must be >= 1.")
    if not 0 < args.min_finite_me_coverage <= 1:
        raise ValueError("--min-finite-me-coverage must be in (0, 1].")
    if args.local_radius <= 0:
        raise ValueError("--local-radius must be positive.")
    if args.composition_radius <= 0:
        raise ValueError("--composition-radius must be positive.")
    if args.randomizations < 999:
        raise ValueError("--randomizations should be at least 999.")
    if args.bootstraps < 500:
        raise ValueError("--bootstraps should be at least 500.")
    if args.matches_per_case < 1:
        raise ValueError("--matches-per-case must be >= 1.")
    if not 1 <= args.min_matches_per_case <= args.matches_per_case:
        raise ValueError(
            "--min-matches-per-case must be between 1 and "
            "--matches-per-case."
        )
    for option_name in (
        "pl_caliper",
        "log_depth_caliper",
        "log_support_caliper",
        "gc_caliper",
        "target_caliper",
    ):
        if getattr(args, option_name) <= 0:
            raise ValueError(f"--{option_name.replace('_', '-')} must be > 0.")
    if args.continuous_block_bp <= 0:
        raise ValueError("--continuous-block-bp must be positive.")

    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}. "
            "Use a new directory to avoid mixing analyses."
        )
    for name in ("qc", "data", "results", "figures", "logs", "work"):
        (output_dir / name).mkdir(parents=True, exist_ok=True)

    configure_logging(output_dir)
    LOG.info("Starting ME-node enrichment version %s.", VERSION)

    contig_lengths = read_fai(args.fai)
    sequences = read_fasta(args.fasta)
    reference = SequenceReference(sequences, contig_lengths)
    annotation = AnnotationIndex(
        args.gff,
        args.operons,
        args.upstream_bp,
    )

    me_inputs: list[tuple[str, str, Path]] = []
    mod_labels: set[str] = set()
    for mod_type, sample, raw_directory in args.me:
        directory = Path(raw_directory).expanduser().resolve()
        require_directory(directory, f"ME directory for {mod_type}/{sample}")
        require_file(directory / "windows.bedgraph", "windows.bedgraph")
        require_file(directory / "regions.bed", "regions.bed")
        me_inputs.append((mod_type, sample, directory))
        mod_labels.add(mod_type)

    provenance_qc = collect_upstream_provenance(me_inputs)
    provenance_qc.to_csv(
        output_dir / "qc" / "upstream_provenance.tsv",
        sep="\t",
        index=False,
    )

    cases_by_mod, treatments, metadata_qc = read_metadata_cases(
        args.metadata,
        mod_labels,
        contig_lengths,
    )
    treatment_to_index = {
        treatment: index
        for index, treatment in enumerate(treatments)
    }
    LOG.info(
        "Parsed %d exact unique metadata case sites: %s.",
        sum(len(group) for group in cases_by_mod.values()),
        ", ".join(
            f"{mod}={len(cases_by_mod.get(mod, {}))}"
            for mod in sorted(mod_labels)
        ),
    )

    unique_case_rows = []
    for mod_type, case_map in cases_by_mod.items():
        for key, case in case_map.items():
            unique_case_rows.append(
                {
                    "mod_type": mod_type,
                    "chrom": key.chrom,
                    "start": key.start,
                    "end": key.end,
                    "strand": key.strand,
                    "treatments": ",".join(sorted(case.treatments)),
                    "treatment_count": len(case.treatments),
                    "direction_from_bracketed_PL": case.direction,
                    "max_abs_treatment_minus_control_PL": (
                        case.max_abs_delta
                    ),
                    "control_PL_by_treatment": json.dumps(
                        {
                            treatment: values
                            for treatment, values
                            in case.control_pl.items()
                        },
                        sort_keys=True,
                    ),
                    "treatment_PL_by_treatment": json.dumps(
                        {
                            treatment: values
                            for treatment, values
                            in case.treatment_pl.items()
                        },
                        sort_keys=True,
                    ),
                    "reported_delta_by_treatment": json.dumps(
                        {
                            treatment: values
                            for treatment, values
                            in case.reported_deltas.items()
                        },
                        sort_keys=True,
                    ),
                    "source_regions": " || ".join(
                        sorted(case.source_regions)
                    ),
                }
            )
    pd.DataFrame(unique_case_rows).to_csv(
        output_dir / "data" / "unique_cases.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    dmr_paths = discover_dmr_files(
        args.dmr_root,
        mod_labels,
        treatments,
    )

    required_dmr_pairs = {
        (mod_type, treatment)
        for mod_type, case_map in cases_by_mod.items()
        for case in case_map.values()
        for treatment in case.treatments
    }
    missing_required_pairs = sorted(
        required_dmr_pairs - set(dmr_paths)
    )
    if missing_required_pairs:
        raise RuntimeError(
            "A metadata case treatment has no matching DMR file: "
            + ", ".join(
                f"{mod}/{treatment}"
                for mod, treatment in missing_required_pairs
            )
        )

    required_dmr_columns = {
        "chrom", "start", "end", "name", "strand",
        "control_total", "treatment_total",
        "control_pct_modified", "treatment_pct_modified",
    }
    for pair, path in dmr_paths.items():
        header = set(dmr_header(path))
        missing_columns = required_dmr_columns - header
        if missing_columns:
            raise ValueError(
                f"{path}: missing required DMR columns "
                f"{sorted(missing_columns)}."
            )

    # Input manifest.
    manifest_rows = []
    core_inputs = [
        ("metadata", "", "", args.metadata),
        ("fasta", "", "", args.fasta),
        ("fai", "", "", args.fai),
        ("gff", "", "", args.gff),
        ("operons", "", "", args.operons),
    ]
    for kind, mod_type, label, path in core_inputs:
        require_file(path, kind)
        manifest_rows.append(
            {
                "type": kind,
                "mod_type": mod_type,
                "label": label,
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
            }
        )
    for (mod_type, treatment), path in sorted(dmr_paths.items()):
        manifest_rows.append(
            {
                "type": "dmr",
                "mod_type": mod_type,
                "label": treatment,
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
            }
        )
    for mod_type, sample, directory in me_inputs:
        for filename in ("windows.bedgraph", "regions.bed"):
            path = directory / filename
            manifest_rows.append(
                {
                    "type": filename,
                    "mod_type": mod_type,
                    "label": sample,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                }
            )
        provenance_candidates = [
            directory / "provenance.json",
            directory / "modkit_entropy.log",
        ]
        for path in provenance_candidates:
            if path.is_file():
                manifest_rows.append(
                    {
                        "type": "provenance",
                        "mod_type": mod_type,
                        "label": sample,
                        "path": str(path),
                        "size_bytes": path.stat().st_size,
                    }
                )

    pd.DataFrame(manifest_rows).to_csv(
        output_dir / "run_manifest.tsv",
        sep="\t",
        index=False,
    )
    metadata_qc.to_csv(
        output_dir / "qc" / "metadata_PL_delta_consistency.tsv",
        sep="\t",
        index=False,
    )

    config = vars(args).copy()
    config["metadata"] = str(args.metadata)
    config["fasta"] = str(args.fasta)
    config["fai"] = str(args.fai)
    config["gff"] = str(args.gff)
    config["operons"] = str(args.operons)
    config["dmr_root"] = str(args.dmr_root)
    config["output_dir"] = str(args.output_dir)
    config["me"] = [
        [mod, sample, str(directory)]
        for mod, sample, directory in me_inputs
    ]
    config["version"] = VERSION
    (output_dir / "analysis_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n"
    )

    pd.DataFrame(
        [
            {
                "question": (
                    "Does stress change ME more at responsive than "
                    "background loci?"
                ),
                "status": "not_testable_from_control_only_ME",
                "required_additional_input": (
                    "control and treatment windows.bedgraph per biological "
                    "replicate"
                ),
            },
            {
                "question": (
                    "Do multi-position patterns show coordination beyond "
                    "their marginal per-position PL values?"
                ),
                "status": "not_testable_from_entropy_summaries",
                "required_additional_input": (
                    "read-level pattern counts or modkit extract output"
                ),
            },
            {
                "question": (
                    "Is the association reproducible among independent "
                    "biological controls?"
                ),
                "status": (
                    "descriptive_only_when_multiple_control_tracks_are_given"
                ),
                "required_additional_input": (
                    "separate control ME tracks with biological replicate IDs"
                ),
            },
        ]
    ).to_csv(
        output_dir / "results" / "questions_not_tested.tsv",
        sep="\t",
        index=False,
    )

    # Parse genuine Modkit windows and region QC.
    tracks: list[WindowTrack] = []
    me_qc_rows: list[dict[str, object]] = []
    region_qc_rows: list[dict[str, object]] = []
    for mod_type, sample, directory in me_inputs:
        sample_tracks, rows = parse_windows_bedgraph(
            mod_type,
            sample,
            directory / "windows.bedgraph",
            contig_lengths,
        )
        tracks.extend(sample_tracks)
        me_qc_rows.extend(rows)
        region_qc_rows.extend(
            parse_regions_qc(
                mod_type,
                sample,
                directory / "regions.bed",
                contig_lengths,
            )
        )

    me_qc = pd.DataFrame(me_qc_rows)
    region_qc = pd.DataFrame(region_qc_rows)
    me_qc.to_csv(
        output_dir / "qc" / "windows_bedgraph_qc.tsv",
        sep="\t",
        index=False,
    )
    region_qc.to_csv(
        output_dir / "qc" / "regions_bed_qc.tsv",
        sep="\t",
        index=False,
    )

    low_window_coverage = me_qc[
        me_qc["finite_coverage_fraction"]
        < args.min_finite_me_coverage
    ]
    low_region_coverage = region_qc[
        region_qc["finite_entropy_coverage_fraction"]
        < args.min_finite_me_coverage
    ]

    # Modkit entropy windows are only emitted where informative modified-base
    # positions and support exist. Genome-wide finite coverage is therefore a
    # warning, not a valid reason to abort a site-level analysis. Formal tests
    # still require finite ME for every retained case and matched control.
    if not low_window_coverage.empty:
        LOG.warning(
            "%d contig/strand ME-window tracks have finite coverage below "
            "the %.1f%% QC warning threshold. Candidate-level availability "
            "and attrition will determine inferential validity. Minimum "
            "observed window coverage was %.1f%%. See "
            "qc/windows_bedgraph_qc.tsv.",
            len(low_window_coverage),
            100.0 * args.min_finite_me_coverage,
            100.0 * float(
                low_window_coverage["finite_coverage_fraction"].min()
            ),
        )

    if not low_region_coverage.empty:
        LOG.warning(
            "%d regions.bed summaries have finite coverage below the %.1f%% "
            "QC warning threshold. regions.bed is QC-only. Minimum observed "
            "regional coverage was %.1f%%. See qc/regions_bed_qc.tsv.",
            len(low_region_coverage),
            100.0 * args.min_finite_me_coverage,
            100.0 * float(
                low_region_coverage[
                    "finite_entropy_coverage_fraction"
                ].min()
            ),
        )

    me_index = MEIndex(tracks)

    # A dry run validates all file schemas but deliberately avoids streaming
    # millions of DMR rows.
    if args.dry_run:
        LOG.info("Dry run completed successfully.")
        write_summary(
            output_dir,
            pd.DataFrame(),
            pd.DataFrame(),
            me_qc,
            metadata_qc,
            me_inputs,
            provenance_qc,
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )
        zip_path = zip_results(output_dir)
        LOG.info("Dry-run report: %s", zip_path)
        return 0

    database_path = output_dir / "work" / "callable_sites.sqlite"
    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")

    for mod_type in mod_labels:
        create_site_table(
            connection,
            mod_type,
            len(treatments),
        )

    dmr_qc_rows = []
    for mod_type in sorted(mod_labels):
        for treatment in treatments:
            path = dmr_paths.get((mod_type, treatment))
            if path is None:
                continue
            LOG.info(
                "Streaming callable sites for %s/%s.",
                mod_type,
                treatment,
            )
            dmr = DMRFile(
                mod_type=mod_type,
                treatment=treatment,
                path=path,
                bit=1 << treatment_to_index[treatment],
                treatment_index=treatment_to_index[treatment],
            )
            dmr_qc_rows.append(
                stream_dmr_into_sqlite(
                    connection,
                    dmr,
                    len(treatments),
                )
            )
    dmr_qc = pd.DataFrame(dmr_qc_rows)
    dmr_qc.to_csv(
        output_dir / "qc" / "dmr_input_qc.tsv",
        sep="\t",
        index=False,
    )

    # A modification label should map consistently to the same canonical base
    # across all of its supplied DMR comparisons. This does not verify the
    # Dorado modification code, but it catches incompatible DMR inputs.
    for mod_type, subset in dmr_qc.groupby("mod_type"):
        bases = {
            str(value)
            for value in subset["dominant_name_base"]
            if str(value) in DNA_BASES
        }
        if len(bases) > 1:
            raise RuntimeError(
                f"DMR files labelled {mod_type!r} contain inconsistent "
                f"dominant canonical bases: {sorted(bases)}."
            )

    # Case callability and coordinate validation.
    case_qc_rows = []
    case_records: dict[
        tuple[str, SiteKey], SiteRecord
    ] = {}
    coordinate_rows = []
    metadata_dmr_pl_rows = []

    all_case_keys_by_mod = {
        mod_type: set(cases)
        for mod_type, cases in cases_by_mod.items()
    }

    for mod_type, cases in cases_by_mod.items():
        for key, case in cases.items():
            record = fetch_site_record(
                connection,
                mod_type,
                key,
                len(treatments),
            )
            required = case_mask(case, treatment_to_index)
            callable_case = (
                record is not None
                and (record.mask & required) == required
            )
            if record is not None:
                case_records[(mod_type, key)] = record
                observed_base = reference.reference_base(key)
                coordinate_rows.append(
                    {
                        "mod_type": mod_type,
                        "chrom": key.chrom,
                        "start": key.start,
                        "end": key.end,
                        "strand": key.strand,
                        "dmr_name": record.name,
                        "strand_oriented_reference_base": observed_base,
                        "one_base_interval": key.end - key.start == 1,
                        "base_matches_DMR_name": (
                            observed_base == record.name
                            if observed_base and record.name in DNA_BASES
                            else np.nan
                        ),
                    }
                )
                for treatment in sorted(case.treatments):
                    treatment_index = treatment_to_index.get(treatment)
                    if treatment_index is None:
                        continue
                    metadata_control_values = case.control_pl.get(
                        treatment,
                        [],
                    )
                    metadata_treatment_values = case.treatment_pl.get(
                        treatment,
                        [],
                    )
                    metadata_control = (
                        float(np.median(metadata_control_values))
                        if metadata_control_values
                        else float("nan")
                    )
                    metadata_treated = (
                        float(np.median(metadata_treatment_values))
                        if metadata_treatment_values
                        else float("nan")
                    )
                    dmr_control = record.pls[treatment_index]
                    dmr_treated = record.post_pls[treatment_index]
                    metadata_dmr_pl_rows.append(
                        {
                            "mod_type": mod_type,
                            "treatment": treatment,
                            "chrom": key.chrom,
                            "start": key.start,
                            "end": key.end,
                            "strand": key.strand,
                            "metadata_control_pl": metadata_control,
                            "dmr_control_pl": dmr_control,
                            "abs_control_difference": (
                                abs(metadata_control - dmr_control)
                                if (
                                    math.isfinite(metadata_control)
                                    and dmr_control is not None
                                )
                                else float("nan")
                            ),
                            "metadata_treatment_pl": metadata_treated,
                            "dmr_treatment_pl": dmr_treated,
                            "abs_treatment_difference": (
                                abs(metadata_treated - dmr_treated)
                                if (
                                    math.isfinite(metadata_treated)
                                    and dmr_treated is not None
                                )
                                else float("nan")
                            ),
                        }
                    )
            case_qc_rows.append(
                {
                    "mod_type": mod_type,
                    "chrom": key.chrom,
                    "start": key.start,
                    "end": key.end,
                    "strand": key.strand,
                    "treatments": ",".join(sorted(case.treatments)),
                    "required_mask": required,
                    "unique_case": 1,
                    "callable": int(callable_case),
                    "scored": 0,
                    "analysed": 0,
                    "direction": case.direction,
                    "max_abs_delta": case.max_abs_delta,
                    "context": annotation.classify(key),
                }
            )

    case_qc = pd.DataFrame(case_qc_rows)
    coordinate_qc = pd.DataFrame(coordinate_rows)
    coordinate_qc.to_csv(
        output_dir / "qc" / "coordinate_and_base_validation.tsv",
        sep="\t",
        index=False,
    )

    metadata_dmr_pl_qc = pd.DataFrame(metadata_dmr_pl_rows)
    metadata_dmr_pl_qc.to_csv(
        output_dir / "qc" / "metadata_vs_DMR_PL.tsv",
        sep="\t",
        index=False,
    )
    if not metadata_dmr_pl_qc.empty:
        differences = pd.concat(
            [
                metadata_dmr_pl_qc["abs_control_difference"],
                metadata_dmr_pl_qc["abs_treatment_difference"],
            ],
            ignore_index=True,
        ).dropna()
        if not differences.empty:
            maximum_difference = float(differences.max())
            if maximum_difference > 0.01:
                raise RuntimeError(
                    "Metadata bracketed PL and DMR PL differ by more than "
                    f"0.01 (maximum={maximum_difference:.6g}). See "
                    "qc/metadata_vs_DMR_PL.tsv."
                )
            if maximum_difference > 1e-5:
                LOG.warning(
                    "Metadata and DMR PL differ slightly "
                    "(maximum absolute difference %.6g).",
                    maximum_difference,
                )

    if not coordinate_qc.empty:
        if not bool(coordinate_qc["one_base_interval"].all()):
            raise RuntimeError(
                "At least one supplied case is not a one-base, zero-based "
                "half-open interval. See qc/coordinate_and_base_validation.tsv."
            )
        comparable = coordinate_qc[
            coordinate_qc["base_matches_DMR_name"].notna()
        ]
        if not comparable.empty:
            mismatch_fraction = float(
                1.0
                - comparable["base_matches_DMR_name"].astype(bool).mean()
            )
            if mismatch_fraction > 0.05:
                raise RuntimeError(
                    f"{mismatch_fraction:.1%} of case coordinates disagree "
                    "with the strand-oriented reference base and DMR name. "
                    "This suggests a coordinate or strand convention mismatch. "
                    "See qc/coordinate_and_base_validation.tsv."
                )
            if mismatch_fraction > 0:
                LOG.warning(
                    "%.2f%% of case coordinates failed reference-base "
                    "validation.",
                    100.0 * mismatch_fraction,
                )

    # Build unbiased, final-mask reservoirs before scoring.
    reservoirs_by_mod = {}
    reservoir_qc_frames = []
    for mod_type in sorted(mod_labels):
        LOG.info("Building random eligible-site reservoir for %s.", mod_type)
        reservoirs, reservoir_qc = build_global_reservoirs(
            connection,
            mod_type,
            len(treatments),
            all_case_keys_by_mod.get(mod_type, set()),
            args.reservoir_per_stratum,
            args.seed,
        )
        reservoirs_by_mod[mod_type] = reservoirs
        reservoir_qc_frames.append(reservoir_qc)

    reservoir_qc = pd.concat(
        reservoir_qc_frames,
        ignore_index=True,
    ) if reservoir_qc_frames else pd.DataFrame()
    reservoir_qc.to_csv(
        output_dir / "qc" / "global_background_sampling.tsv",
        sep="\t",
        index=False,
    )

    # Score all callable cases and all reservoir controls at all requested
    # support thresholds. Local controls are scored lazily and cached.
    support_thresholds = sorted(
        {
            args.min_window_support,
            *args.support_sensitivity,
        }
    )
    score_cache: dict[
        tuple[str, SiteKey, int], Optional[SiteScore]
    ] = {}

    unique_reservoir_records: dict[
        tuple[str, SiteKey], SiteRecord
    ] = {}
    for mod_type, reservoirs in reservoirs_by_mod.items():
        for reservoir in reservoirs.values():
            for record in reservoir.items:
                unique_reservoir_records[(mod_type, record.key)] = record

    total_to_score = (
        len(case_records)
        + len(unique_reservoir_records)
    ) * len(support_thresholds)
    LOG.info(
        "Scoring %d site/threshold combinations from true windows.bedgraph.",
        total_to_score,
    )

    for (mod_type, key), _ in {
        **case_records,
        **unique_reservoir_records,
    }.items():
        for threshold in support_thresholds:
            score_cache[(mod_type, key, threshold)] = me_index.score(
                mod_type,
                key,
                threshold,
            )

    # Build primary case-specific risk sets.
    analyses_by_mod: dict[str, list[CaseAnalysis]] = defaultdict(list)
    case_level_rows = []
    control_pool_rows = []
    matched_set_rows = []
    lazy_local_seen: dict[
        tuple[str, SiteKey], SiteRecord
    ] = {}

    for row_index, qc_row in case_qc.iterrows():
        mod_type = str(qc_row["mod_type"])
        key = SiteKey(
            str(qc_row["chrom"]),
            int(qc_row["start"]),
            int(qc_row["end"]),
            str(qc_row["strand"]),
        )
        case = cases_by_mod[mod_type][key]
        record = case_records.get((mod_type, key))
        score = score_cache.get(
            (mod_type, key, args.min_window_support)
        )
        if record is None or score is None:
            continue
        case_qc.loc[row_index, "scored"] = 1

        local_records = query_local_controls(
            connection,
            mod_type,
            case,
            len(treatments),
            treatment_to_index,
            args.local_radius,
            all_case_keys_by_mod[mod_type],
        )
        for local_record in local_records:
            lazy_local_seen[(mod_type, local_record.key)] = local_record
            cache_key = (
                mod_type,
                local_record.key,
                args.min_window_support,
            )
            if cache_key not in score_cache:
                score_cache[cache_key] = me_index.score(
                    mod_type,
                    local_record.key,
                    args.min_window_support,
                )

        global_records = eligible_global_controls(
            reservoirs_by_mod[mod_type],
            case,
            treatment_to_index,
        )

        domain_endpoint = f"domain_{args.domain_radius}"
        analysis = make_case_analysis(
            case=case,
            case_record=record,
            case_score=score,
            local_records=local_records,
            global_records=global_records,
            score_cache=score_cache,
            reference=reference,
            composition_radius=args.composition_radius,
            min_support=args.min_window_support,
            treatment_to_index=treatment_to_index,
            min_local_controls=args.min_local_controls,
            min_global_controls=args.min_global_controls,
            matches_per_case=args.matches_per_case,
            min_matches_per_case=args.min_matches_per_case,
            pl_caliper=args.pl_caliper,
            log_depth_caliper=args.log_depth_caliper,
            log_support_caliper=args.log_support_caliper,
            gc_caliper=args.gc_caliper,
            target_caliper=args.target_caliper,
            domain_endpoint=domain_endpoint,
        )
        if analysis is None:
            continue

        analyses_by_mod[mod_type].append(analysis)
        case_qc.loc[row_index, "analysed"] = 1

        for endpoint_name, control_keys, control_values, case_value in (
            (
                "local",
                analysis.local_control_keys,
                analysis.local_control_values,
                score.site_me,
            ),
            (
                "global",
                analysis.global_control_keys,
                analysis.global_control_values,
                score.endpoint(domain_endpoint),
            ),
        ):
            for control_rank, (control_key, control_value) in enumerate(
                zip(control_keys, control_values),
                start=1,
            ):
                matched_set_rows.append(
                    {
                        "mod_type": mod_type,
                        "endpoint": endpoint_name,
                        "case_chrom": key.chrom,
                        "case_start": key.start,
                        "case_end": key.end,
                        "case_strand": key.strand,
                        "control_rank": control_rank,
                        "control_chrom": control_key.chrom,
                        "control_start": control_key.start,
                        "control_end": control_key.end,
                        "control_strand": control_key.strand,
                        "case_ME": case_value,
                        "control_ME": float(control_value),
                        "case_minus_control_ME": (
                            case_value - float(control_value)
                        ),
                        "matching_used_outcome": False,
                    }
                )

        case_level_rows.append(
            {
                "mod_type": mod_type,
                "chrom": key.chrom,
                "start": key.start,
                "end": key.end,
                "strand": key.strand,
                "strand_mode": score.strand_mode,
                "treatments": ",".join(sorted(case.treatments)),
                "direction": case.direction,
                "max_abs_delta": case.max_abs_delta,
                "context": annotation.classify(key),
                "site_me": score.site_me,
                "site_q75": score.site_q75,
                "site_support": score.site_support,
                "overlap_windows": score.overlap_windows,
                "domain_me": score.endpoint(domain_endpoint),
                "context_500_me": score.context_500_excluding_site,
                "local_prominence_500": score.local_prominence_500,
                "local_percentile_500": score.local_percentile_500,
                "local_peak90": (
                    float(score.local_percentile_500 >= 0.90)
                    if math.isfinite(score.local_percentile_500)
                    else float("nan")
                ),
                "local_peak95": (
                    float(score.local_percentile_500 >= 0.95)
                    if math.isfinite(score.local_percentile_500)
                    else float("nan")
                ),
                "local_window_prominence": score.local_window_prominence,
                "local_effect_vs_callable_neighbours": analysis.local_effect,
                "local_percentile": analysis.local_percentile,
                "global_effect_vs_random_eligible_sites": analysis.global_effect,
                "global_percentile": analysis.global_percentile,
                "local_control_count": analysis.local_controls,
                "local_control_ESS": analysis.local_ess,
                "global_control_count": analysis.global_controls,
                "global_control_ESS": analysis.global_ess,
                "local_within_case_max_standardized_mismatch": (
                    analysis.local_within_case_max_standardized_mismatch
                ),
                "global_within_case_max_standardized_mismatch": (
                    analysis.global_within_case_max_standardized_mismatch
                ),
                "case_control_pl": analysis.case_control_pl,
                "case_control_depth": analysis.case_control_depth,
                "case_me_support": analysis.case_me_support,
                "case_gc_fraction": analysis.case_gc_fraction,
                "case_target_base_fraction": (
                    analysis.case_target_base_fraction
                ),
                "local_expected_control_pl": (
                    analysis.local_expected_control_pl
                ),
                "local_expected_control_depth": (
                    analysis.local_expected_control_depth
                ),
                "local_expected_me_support": (
                    analysis.local_expected_me_support
                ),
                "local_expected_gc_fraction": (
                    analysis.local_expected_gc_fraction
                ),
                "local_expected_target_base_fraction": (
                    analysis.local_expected_target_base_fraction
                ),
                "global_expected_control_pl": (
                    analysis.global_expected_control_pl
                ),
                "global_expected_control_depth": (
                    analysis.global_expected_control_depth
                ),
                "global_expected_me_support": (
                    analysis.global_expected_me_support
                ),
                "global_expected_gc_fraction": (
                    analysis.global_expected_gc_fraction
                ),
                "global_expected_target_base_fraction": (
                    analysis.global_expected_target_base_fraction
                ),
            }
        )
        control_pool_rows.append(
            {
                "mod_type": mod_type,
                "chrom": key.chrom,
                "start": key.start,
                "strand": key.strand,
                "local_controls": analysis.local_controls,
                "local_ESS": analysis.local_ess,
                "global_controls": analysis.global_controls,
                "global_ESS": analysis.global_ess,
                "local_within_case_max_standardized_mismatch": (
                    analysis.local_within_case_max_standardized_mismatch
                ),
                "global_within_case_max_standardized_mismatch": (
                    analysis.global_within_case_max_standardized_mismatch
                ),
                "case_control_pl": analysis.case_control_pl,
                "case_control_depth": analysis.case_control_depth,
                "case_me_support": analysis.case_me_support,
                "case_gc_fraction": analysis.case_gc_fraction,
                "case_target_base_fraction": (
                    analysis.case_target_base_fraction
                ),
                "local_expected_control_pl": (
                    analysis.local_expected_control_pl
                ),
                "local_expected_control_depth": (
                    analysis.local_expected_control_depth
                ),
                "local_expected_me_support": (
                    analysis.local_expected_me_support
                ),
                "local_expected_gc_fraction": (
                    analysis.local_expected_gc_fraction
                ),
                "local_expected_target_base_fraction": (
                    analysis.local_expected_target_base_fraction
                ),
                "global_expected_control_pl": (
                    analysis.global_expected_control_pl
                ),
                "global_expected_control_depth": (
                    analysis.global_expected_control_depth
                ),
                "global_expected_me_support": (
                    analysis.global_expected_me_support
                ),
                "global_expected_gc_fraction": (
                    analysis.global_expected_gc_fraction
                ),
                "global_expected_target_base_fraction": (
                    analysis.global_expected_target_base_fraction
                ),
            }
        )

    case_qc.to_csv(
        output_dir / "qc" / "case_attrition.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(control_pool_rows).to_csv(
        output_dir / "qc" / "case_specific_control_pools.tsv",
        sep="\t",
        index=False,
    )

    balance_qc = make_balance_table(analyses_by_mod)
    if not balance_qc.empty:
        balance_qc["stage"] = "before_outcome_blind_balance_trim"
    balance_qc.to_csv(
        output_dir / "qc" / "covariate_balance.tsv",
        sep="\t",
        index=False,
    )

    case_level = pd.DataFrame(case_level_rows)
    case_level.to_csv(
        output_dir / "data" / "case_level_results.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(matched_set_rows).to_csv(
        output_dir / "data" / "matched_sets.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Primary local and global tests.
    primary_rows = []
    for mod_type in sorted(mod_labels):
        analyses = analyses_by_mod.get(mod_type, [])
        if len(analyses) < 10:
            for endpoint in ("local", "global"):
                primary_rows.append(
                    {
                        "endpoint": endpoint,
                        "site_cases": len(analyses),
                        "event_clusters": len(
                            spatial_clusters(analyses, args.cluster_bp)
                        ) if analyses else 0,
                        "cluster_bp": args.cluster_bp,
                        "minimum_cluster_rule_passed": False,
                        "inference_validity": (
                            "underpowered_fewer_than_10_valid_cases"
                        ),
                        "observed_mean_effect": float("nan"),
                        "bootstrap_ci_low": float("nan"),
                        "bootstrap_ci_high": float("nan"),
                        "genomic_randomization_p_high": float("nan"),
                        "genomic_randomization_p_low": float("nan"),
                        "genomic_randomization_p_two_sided": float("nan"),
                        "tail90_randomization_p_high": float("nan"),
                        "tail95_randomization_p_high": float("nan"),
                        "mod_type": mod_type,
                        "analysis": (
                            "local_site_specificity"
                            if endpoint == "local"
                            else "global_domain_enrichment"
                        ),
                        "support_threshold": args.min_window_support,
                        "local_radius_bp": args.local_radius,
                        "domain_radius_bp": args.domain_radius,
                        "composition_radius_bp": args.composition_radius,
                    }
                )
            continue
        for endpoint in ("local", "global"):
            result = randomization_test(
                analyses,
                endpoint,
                randomizations=args.randomizations,
                bootstraps=args.bootstraps,
                cluster_bp=args.cluster_bp,
                seed=stable_seed(
                    args.seed,
                    "primary",
                    mod_type,
                    endpoint,
                ),
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "analysis": (
                        "local_site_specificity"
                        if endpoint == "local"
                        else "global_domain_enrichment"
                    ),
                    "support_threshold": args.min_window_support,
                    "local_radius_bp": args.local_radius,
                    "domain_radius_bp": args.domain_radius,
                    "composition_radius_bp": args.composition_radius,
                }
            )
            primary_rows.append(result)

    primary = pd.DataFrame(primary_rows)
    if not primary.empty:
        # Balance is evaluated inside every randomization test, including
        # every sensitivity configuration.
        primary["holm_primary_p_high"] = holm_adjust(
            primary["genomic_randomization_p_high"].tolist()
        )
        primary["holm_primary_p_low"] = holm_adjust(
            primary["genomic_randomization_p_low"].tolist()
        )
        primary["holm_primary_p_two_sided"] = holm_adjust(
            primary["genomic_randomization_p_two_sided"].tolist()
        )
        # Backward-friendly alias: the prespecified biological hypothesis is
        # enrichment toward higher control ME.
        primary["holm_primary_p"] = primary["holm_primary_p_high"]
        primary["interpretation"] = [
            (
                interpretation_label(
                    float(row["observed_mean_effect"]),
                    float(row["holm_primary_p_high"]),
                    float(row["holm_primary_p_low"]),
                    float(row["holm_primary_p_two_sided"]),
                )
                if row["inference_validity"]
                == "valid_matched_set_genomic_test"
                else str(row["inference_validity"])
            )
            for _, row in primary.iterrows()
        ]
    primary.to_csv(
        output_dir / "results" / "primary_tests.tsv",
        sep="\t",
        index=False,
    )
    if not primary.empty:
        posttrim_columns = [
            column
            for column in (
                "mod_type",
                "endpoint",
                "site_cases_before_balance_trim",
                "balance_trimmed_cases",
                "balance_trim_fraction",
                "balance_target_achieved",
                "smd_control_pl",
                "smd_log_control_depth",
                "smd_log_me_support",
                "smd_gc_fraction",
                "smd_target_base_fraction",
                "max_abs_covariate_smd",
                "inference_validity",
            )
            if column in primary.columns
        ]
        primary[posttrim_columns].to_csv(
            output_dir / "qc" / "primary_posttrim_balance.tsv",
            sep="\t",
            index=False,
        )

    # Direction and treatment analyses.
    direction_rows = []
    treatment_rows = []
    for mod_type, analyses in analyses_by_mod.items():
        for direction in ("gain", "loss"):
            subset = [
                analysis
                for analysis in analyses
                if analysis.case.direction == direction
            ]
            if len(subset) < 10:
                continue
            for endpoint in ("local", "global"):
                result = randomization_test(
                    subset,
                    endpoint,
                    randomizations=args.randomizations,
                    bootstraps=args.bootstraps,
                    cluster_bp=args.cluster_bp,
                    seed=stable_seed(
                        args.seed,
                        "direction",
                        mod_type,
                        direction,
                        endpoint,
                    ),
                )
                result.update(
                    {
                        "mod_type": mod_type,
                        "direction": direction,
                    }
                )
                direction_rows.append(result)

        gain = [
            analysis
            for analysis in analyses
            if analysis.case.direction == "gain"
        ]
        loss = [
            analysis
            for analysis in analyses
            if analysis.case.direction == "loss"
        ]
        for endpoint in ("local", "global"):
            interaction = direction_difference_test(
                gain,
                loss,
                endpoint,
                randomizations=args.randomizations,
                bootstraps=args.bootstraps,
                cluster_bp=args.cluster_bp,
                seed=stable_seed(
                    args.seed,
                    "direction_interaction",
                    mod_type,
                    endpoint,
                ),
            )
            interaction.update(
                {
                    "mod_type": mod_type,
                    "direction": "gain_minus_loss",
                }
            )
            direction_rows.append(interaction)

        for treatment in treatments:
            subset = [
                analysis
                for analysis in analyses
                if treatment in analysis.case.treatments
            ]
            if len(subset) < 10:
                continue
            for endpoint in ("local", "global"):
                result = randomization_test(
                    subset,
                    endpoint,
                    randomizations=args.randomizations,
                    bootstraps=args.bootstraps,
                    cluster_bp=args.cluster_bp,
                    seed=stable_seed(
                        args.seed,
                        "treatment",
                        mod_type,
                        treatment,
                        endpoint,
                    ),
                )
                result.update(
                    {
                        "mod_type": mod_type,
                        "treatment": treatment,
                    }
                )
                treatment_rows.append(result)

    direction_frame = pd.DataFrame(direction_rows)
    if not direction_frame.empty:
        if "genomic_randomization_p_high" in direction_frame:
            direction_frame["bh_high_p"] = bh_adjust(
                direction_frame[
                    "genomic_randomization_p_high"
                ].fillna(1.0).tolist()
            )
            direction_frame["bh_low_p"] = bh_adjust(
                direction_frame[
                    "genomic_randomization_p_low"
                ].fillna(1.0).tolist()
            )
            direction_frame["bh_two_sided_p"] = bh_adjust(
                direction_frame[
                    "genomic_randomization_p_two_sided"
                ].fillna(1.0).tolist()
            )
            direction_frame["bh_tail90_p"] = bh_adjust(
                direction_frame[
                    "tail90_randomization_p_high"
                ].fillna(1.0).tolist()
            )
            direction_frame["bh_tail95_p"] = bh_adjust(
                direction_frame[
                    "tail95_randomization_p_high"
                ].fillna(1.0).tolist()
            )
            direction_frame["bh_p"] = direction_frame["bh_high_p"]
        if "randomization_p_two_sided" in direction_frame:
            direction_frame["bh_interaction_two_sided_p"] = bh_adjust(
                direction_frame[
                    "randomization_p_two_sided"
                ].fillna(1.0).tolist()
            )
    direction_frame.to_csv(
        output_dir / "results" / "direction_tests.tsv",
        sep="\t",
        index=False,
    )

    treatment_frame = pd.DataFrame(treatment_rows)
    if not treatment_frame.empty:
        treatment_frame["bh_high_p"] = bh_adjust(
            treatment_frame[
                "genomic_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        treatment_frame["bh_low_p"] = bh_adjust(
            treatment_frame[
                "genomic_randomization_p_low"
            ].fillna(1.0).tolist()
        )
        treatment_frame["bh_two_sided_p"] = bh_adjust(
            treatment_frame[
                "genomic_randomization_p_two_sided"
            ].fillna(1.0).tolist()
        )
        treatment_frame["bh_tail90_p"] = bh_adjust(
            treatment_frame[
                "tail90_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        treatment_frame["bh_tail95_p"] = bh_adjust(
            treatment_frame[
                "tail95_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        treatment_frame["bh_p"] = treatment_frame["bh_high_p"]
    treatment_frame.to_csv(
        output_dir / "results" / "treatment_tests.tsv",
        sep="\t",
        index=False,
    )

    # Spatial clustering sensitivity on the already valid primary risk sets.
    sensitivity_rows = []
    for mod_type, analyses in analyses_by_mod.items():
        if len(analyses) < 10:
            continue
        for cluster_bp in (0, 100, 200, 500):
            for endpoint in ("local", "global"):
                result = randomization_test(
                    analyses,
                    endpoint,
                    randomizations=max(5000, args.randomizations // 4),
                    bootstraps=max(1000, args.bootstraps // 2),
                    cluster_bp=cluster_bp,
                    seed=stable_seed(
                        args.seed,
                        "cluster_sensitivity",
                        mod_type,
                        endpoint,
                        cluster_bp,
                    ),
                )
                result.update(
                    {
                        "mod_type": mod_type,
                        "sensitivity": "cluster_bp",
                        "sensitivity_value": cluster_bp,
                    }
                )
                sensitivity_rows.append(result)

    # Matching-caliper sensitivity. Outcome values are never used when
    # rebuilding these stricter/looser matched sets.
    for caliper_multiplier in (0.75, 1.25):
        for mod_type in sorted(mod_labels):
            analyses = build_analyses_for_configuration(
                connection=connection,
                mod_type=mod_type,
                cases=cases_by_mod.get(mod_type, {}),
                case_records=case_records,
                all_case_keys=all_case_keys_by_mod.get(mod_type, set()),
                reservoirs=reservoirs_by_mod.get(mod_type, {}),
                score_cache=score_cache,
                me_index=me_index,
                annotation=annotation,
                reference=reference,
                composition_radius=args.composition_radius,
                treatment_count=len(treatments),
                treatment_to_index=treatment_to_index,
                min_support=args.min_window_support,
                local_radius=args.local_radius,
                domain_radius=args.domain_radius,
                min_local_controls=args.min_local_controls,
                min_global_controls=args.min_global_controls,
                matches_per_case=args.matches_per_case,
                min_matches_per_case=args.min_matches_per_case,
                pl_caliper=args.pl_caliper * caliper_multiplier,
                log_depth_caliper=(
                    args.log_depth_caliper * caliper_multiplier
                ),
                log_support_caliper=(
                    args.log_support_caliper * caliper_multiplier
                ),
                gc_caliper=args.gc_caliper * caliper_multiplier,
                target_caliper=(
                    args.target_caliper * caliper_multiplier
                ),
            )
            if len(analyses) < 10:
                continue
            for endpoint in ("local", "global"):
                result = randomization_test(
                    analyses,
                    endpoint,
                    randomizations=max(
                        5000,
                        args.randomizations // 4,
                    ),
                    bootstraps=max(
                        1000,
                        args.bootstraps // 2,
                    ),
                    cluster_bp=args.cluster_bp,
                    seed=stable_seed(
                        args.seed,
                        "caliper_sensitivity",
                        mod_type,
                        endpoint,
                        caliper_multiplier,
                    ),
                )
                result.update(
                    {
                        "mod_type": mod_type,
                        "sensitivity": "caliper_multiplier",
                        "sensitivity_value": caliper_multiplier,
                    }
                )
                sensitivity_rows.append(result)

    # Minimum-support sensitivity. Entropy estimates from very small read
    # counts are unstable, so the conclusion must not depend on one threshold.
    for threshold in support_thresholds:
        if threshold == args.min_window_support:
            continue
        for mod_type in sorted(mod_labels):
            analyses = build_analyses_for_configuration(
                connection=connection,
                mod_type=mod_type,
                cases=cases_by_mod.get(mod_type, {}),
                case_records=case_records,
                all_case_keys=all_case_keys_by_mod.get(mod_type, set()),
                reservoirs=reservoirs_by_mod.get(mod_type, {}),
                score_cache=score_cache,
                me_index=me_index,
                annotation=annotation,
                reference=reference,
                composition_radius=args.composition_radius,
                treatment_count=len(treatments),
                treatment_to_index=treatment_to_index,
                min_support=threshold,
                local_radius=args.local_radius,
                domain_radius=args.domain_radius,
                min_local_controls=args.min_local_controls,
                min_global_controls=args.min_global_controls,
                matches_per_case=args.matches_per_case,
                min_matches_per_case=args.min_matches_per_case,
                pl_caliper=args.pl_caliper,
                log_depth_caliper=args.log_depth_caliper,
                log_support_caliper=args.log_support_caliper,
                gc_caliper=args.gc_caliper,
                target_caliper=args.target_caliper,
            )
            if len(analyses) < 10:
                continue
            for endpoint in ("local", "global"):
                result = randomization_test(
                    analyses,
                    endpoint,
                    randomizations=max(5000, args.randomizations // 4),
                    bootstraps=max(1000, args.bootstraps // 2),
                    cluster_bp=args.cluster_bp,
                    seed=stable_seed(
                        args.seed,
                        "support_sensitivity",
                        mod_type,
                        threshold,
                        endpoint,
                    ),
                )
                result.update(
                    {
                        "mod_type": mod_type,
                        "sensitivity": "minimum_window_support",
                        "sensitivity_value": threshold,
                    }
                )
                sensitivity_rows.append(result)

    # Domain-radius sensitivity for the genome-wide regional question.
    for radius in (250, 1000):
        if radius == args.domain_radius:
            continue
        for mod_type in sorted(mod_labels):
            analyses = build_analyses_for_configuration(
                connection=connection,
                mod_type=mod_type,
                cases=cases_by_mod.get(mod_type, {}),
                case_records=case_records,
                all_case_keys=all_case_keys_by_mod.get(mod_type, set()),
                reservoirs=reservoirs_by_mod.get(mod_type, {}),
                score_cache=score_cache,
                me_index=me_index,
                annotation=annotation,
                reference=reference,
                composition_radius=args.composition_radius,
                treatment_count=len(treatments),
                treatment_to_index=treatment_to_index,
                min_support=args.min_window_support,
                local_radius=args.local_radius,
                domain_radius=radius,
                min_local_controls=args.min_local_controls,
                min_global_controls=args.min_global_controls,
                matches_per_case=args.matches_per_case,
                min_matches_per_case=args.min_matches_per_case,
                pl_caliper=args.pl_caliper,
                log_depth_caliper=args.log_depth_caliper,
                log_support_caliper=args.log_support_caliper,
                gc_caliper=args.gc_caliper,
                target_caliper=args.target_caliper,
            )
            if len(analyses) < 10:
                continue
            result = randomization_test(
                analyses,
                "global",
                randomizations=max(5000, args.randomizations // 4),
                bootstraps=max(1000, args.bootstraps // 2),
                cluster_bp=args.cluster_bp,
                seed=stable_seed(
                    args.seed,
                    "domain_radius_sensitivity",
                    mod_type,
                    radius,
                ),
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "sensitivity": "domain_radius_bp",
                    "sensitivity_value": radius,
                }
            )
            sensitivity_rows.append(result)

    # Local-neighbourhood-radius sensitivity for the site-specific question.
    for radius in (1000, 5000):
        if radius == args.local_radius:
            continue
        for mod_type in sorted(mod_labels):
            analyses = build_analyses_for_configuration(
                connection=connection,
                mod_type=mod_type,
                cases=cases_by_mod.get(mod_type, {}),
                case_records=case_records,
                all_case_keys=all_case_keys_by_mod.get(mod_type, set()),
                reservoirs=reservoirs_by_mod.get(mod_type, {}),
                score_cache=score_cache,
                me_index=me_index,
                annotation=annotation,
                reference=reference,
                composition_radius=args.composition_radius,
                treatment_count=len(treatments),
                treatment_to_index=treatment_to_index,
                min_support=args.min_window_support,
                local_radius=radius,
                domain_radius=args.domain_radius,
                min_local_controls=args.min_local_controls,
                min_global_controls=args.min_global_controls,
                matches_per_case=args.matches_per_case,
                min_matches_per_case=args.min_matches_per_case,
                pl_caliper=args.pl_caliper,
                log_depth_caliper=args.log_depth_caliper,
                log_support_caliper=args.log_support_caliper,
                gc_caliper=args.gc_caliper,
                target_caliper=args.target_caliper,
            )
            if len(analyses) < 10:
                continue
            result = randomization_test(
                analyses,
                "local",
                randomizations=max(5000, args.randomizations // 4),
                bootstraps=max(1000, args.bootstraps // 2),
                cluster_bp=args.cluster_bp,
                seed=stable_seed(
                    args.seed,
                    "local_radius_sensitivity",
                    mod_type,
                    radius,
                ),
            )
            result.update(
                {
                    "mod_type": mod_type,
                    "sensitivity": "local_radius_bp",
                    "sensitivity_value": radius,
                }
            )
            sensitivity_rows.append(result)

    # Optional annotation-matched global sensitivity. Annotation is not used
    # in the primary test, preventing organism-specific annotation choices from
    # defining the result.
    for mod_type in sorted(mod_labels):
        analyses = build_analyses_for_configuration(
            connection=connection,
            mod_type=mod_type,
            cases=cases_by_mod.get(mod_type, {}),
            case_records=case_records,
            all_case_keys=all_case_keys_by_mod.get(mod_type, set()),
            reservoirs=reservoirs_by_mod.get(mod_type, {}),
            score_cache=score_cache,
            me_index=me_index,
            annotation=annotation,
            reference=reference,
            composition_radius=args.composition_radius,
            treatment_count=len(treatments),
            treatment_to_index=treatment_to_index,
            min_support=args.min_window_support,
            local_radius=args.local_radius,
            domain_radius=args.domain_radius,
            min_local_controls=args.min_local_controls,
            min_global_controls=args.min_global_controls,
            matches_per_case=args.matches_per_case,
            min_matches_per_case=args.min_matches_per_case,
            pl_caliper=args.pl_caliper,
            log_depth_caliper=args.log_depth_caliper,
            log_support_caliper=args.log_support_caliper,
            gc_caliper=args.gc_caliper,
            target_caliper=args.target_caliper,
            context_match_global=True,
        )
        if len(analyses) < 10:
            continue
        result = randomization_test(
            analyses,
            "global",
            randomizations=max(5000, args.randomizations // 4),
            bootstraps=max(1000, args.bootstraps // 2),
            cluster_bp=args.cluster_bp,
            seed=stable_seed(
                args.seed,
                "context_matched_sensitivity",
                mod_type,
            ),
        )
        result.update(
            {
                "mod_type": mod_type,
                "sensitivity": "annotation_context_matched",
                "sensitivity_value": "yes",
            }
        )
        sensitivity_rows.append(result)

    sensitivity_frame = pd.DataFrame(sensitivity_rows)
    if not sensitivity_frame.empty:
        sensitivity_frame["bh_high_p"] = bh_adjust(
            sensitivity_frame[
                "genomic_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        sensitivity_frame["bh_low_p"] = bh_adjust(
            sensitivity_frame[
                "genomic_randomization_p_low"
            ].fillna(1.0).tolist()
        )
        sensitivity_frame["bh_two_sided_p"] = bh_adjust(
            sensitivity_frame[
                "genomic_randomization_p_two_sided"
            ].fillna(1.0).tolist()
        )
        sensitivity_frame["bh_tail90_p"] = bh_adjust(
            sensitivity_frame[
                "tail90_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        sensitivity_frame["bh_tail95_p"] = bh_adjust(
            sensitivity_frame[
                "tail95_randomization_p_high"
            ].fillna(1.0).tolist()
        )
        sensitivity_frame["bh_p"] = sensitivity_frame["bh_high_p"]
    sensitivity_frame.to_csv(
        output_dir / "results" / "sensitivity_tests.tsv",
        sep="\t",
        index=False,
    )

    # Continuous response analysis across the random callable risk-set sample.
    # This avoids reducing the biological response to significant/non-significant
    # case status and asks whether baseline ME predicts subsequent |delta PL|
    # after flexible adjustment for baseline PL and technical covariates.
    records_by_mod_for_continuous: dict[
        str,
        dict[SiteKey, SiteRecord],
    ] = defaultdict(dict)
    for (mod_type, key), record in unique_reservoir_records.items():
        records_by_mod_for_continuous[mod_type][key] = record
    for (mod_type, key), record in case_records.items():
        records_by_mod_for_continuous[mod_type][key] = record

    continuous_data = build_continuous_response_table(
        records_by_mod=records_by_mod_for_continuous,
        cases_by_mod=cases_by_mod,
        score_cache=score_cache,
        reference=reference,
        annotation=annotation,
        treatments=treatments,
        treatment_to_index=treatment_to_index,
        min_support=args.min_window_support,
        domain_radius=args.domain_radius,
        composition_radius=args.composition_radius,
        block_bp=args.continuous_block_bp,
    )

    # Attach explicit genome-background references. These distinguish:
    # (i) windows physically overlapping the base;
    # (ii) the non-overlapping +/-500-bp context; and
    # (iii) the callable genome-wide background.
    continuous_data = add_genome_reference_columns(continuous_data)
    union_site_data = build_union_site_frame(continuous_data)

    continuous_data.to_csv(
        output_dir / "data" / "continuous_response_sites.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    union_site_data.to_csv(
        output_dir / "data" / "union_site_risk_set.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Preserve the Version-4 adjusted continuous models and add transparent
    # unadjusted counterparts.
    continuous_adjusted = run_continuous_response_models(
        continuous_data
    )
    continuous_unadjusted = run_unadjusted_continuous_models(
        continuous_data
    )
    continuous_results = pd.concat(
        [continuous_adjusted, continuous_unadjusted],
        ignore_index=True,
        sort=False,
    )
    continuous_results.to_csv(
        output_dir / "results" / "continuous_response_models.tsv",
        sep="\t",
        index=False,
    )

    # Inclusive threshold-positive analyses. Unlike strict matching, these
    # retain every case with valid ME and estimate both adjusted and raw
    # associations across the callable risk set.
    binary_results = run_binary_significance_models(
        continuous_data,
        union_site_data,
    )
    binary_results.to_csv(
        output_dir / "results" / "binary_significance_models.tsv",
        sep="\t",
        index=False,
    )

    case_me_differences = run_case_me_difference_models(
        continuous_data,
        union_site_data,
    )
    case_me_differences.to_csv(
        output_dir / "results" / "case_ME_difference_models.tsv",
        sep="\t",
        index=False,
    )

    direct_peak_contrasts = run_direct_peak_contrasts(
        union_site_data
    )
    direct_peak_contrasts.to_csv(
        output_dir / "results" / "direct_peak_contrasts.tsv",
        sep="\t",
        index=False,
    )

    # Compact machine-readable table of primary and notable exploratory rows.
    key_frames: list[pd.DataFrame] = []
    if not primary.empty:
        frame = primary.copy()
        frame["result_source"] = "primary"
        frame["evidence_level"] = "prespecified_primary"
        key_frames.append(frame)

    if continuous_results is not None and not continuous_results.empty:
        pooled_continuous = continuous_results[
            continuous_results["analysis_scope"]
            == "pooled_treatment_adjusted"
        ].copy()
        if not pooled_continuous.empty:
            pooled_continuous["result_source"] = "continuous_response"
            pooled_continuous["evidence_level"] = (
                "prespecified_secondary_continuous"
            )
            key_frames.append(pooled_continuous)

        significant_continuous = continuous_results[
            (
                continuous_results["analysis_scope"]
                == "treatment_specific"
            )
            & (
                continuous_results[
                    "bh_p_treatment_specific"
                ].fillna(1.0)
                < 0.05
            )
        ].copy()
        if not significant_continuous.empty:
            significant_continuous["result_source"] = (
                "continuous_response"
            )
            significant_continuous["evidence_level"] = (
                "multiplicity_adjusted_exploratory"
            )
            key_frames.append(significant_continuous)

    for source_name, frame in (
        ("direction", direction_frame),
        ("treatment", treatment_frame),
        ("sensitivity", sensitivity_frame),
    ):
        if frame is None or frame.empty:
            continue
        include = pd.Series(False, index=frame.index)

        if "bh_high_p" in frame:
            include |= frame["bh_high_p"].fillna(1.0) < 0.05
        if (
            "bh_low_p" in frame
            and "bh_two_sided_p" in frame
        ):
            include |= (
                (frame["bh_low_p"].fillna(1.0) < 0.05)
                & (
                    frame["bh_two_sided_p"].fillna(1.0)
                    < 0.05
                )
            )
        if "bh_interaction_two_sided_p" in frame:
            include |= (
                frame["bh_interaction_two_sided_p"]
                .fillna(1.0)
                < 0.05
            )

        adjusted = frame.loc[include].copy()
        if not adjusted.empty:
            adjusted["result_source"] = source_name
            adjusted["evidence_level"] = "multiplicity_adjusted_exploratory"
            key_frames.append(adjusted)

        nominal = pd.Series(False, index=frame.index)
        if "genomic_randomization_p_high" in frame:
            nominal |= (
                frame["genomic_randomization_p_high"]
                .fillna(1.0)
                < 0.05
            )
        if "genomic_randomization_p_low" in frame:
            nominal |= (
                frame["genomic_randomization_p_low"]
                .fillna(1.0)
                < 0.05
            )
        if "randomization_p_two_sided" in frame:
            nominal |= (
                frame["randomization_p_two_sided"]
                .fillna(1.0)
                < 0.05
            )
        nominal &= ~include
        nominal_frame = frame.loc[nominal].copy()
        if not nominal_frame.empty:
            nominal_frame["result_source"] = source_name
            nominal_frame["evidence_level"] = "nominal_unadjusted_only"
            key_frames.append(nominal_frame)

    key_results = (
        pd.concat(key_frames, ignore_index=True, sort=False)
        if key_frames
        else pd.DataFrame()
    )
    key_results.to_csv(
        output_dir / "results" / "key_results.tsv",
        sep="\t",
        index=False,
    )

    # Per-control score summaries. With one pooled file this makes the
    # limitation explicit; with repeated --me entries it reveals consistency.
    per_control_rows = []
    for mod_type, analyses in analyses_by_mod.items():
        for sample in me_index.samples_by_mod.get(mod_type, []):
            # A compact descriptive effect based on candidate overlapping-window
            # scores in this control track. No sample-level p-value is invented.
            values = []
            for analysis in analyses:
                track, _ = me_index.track_for(
                    mod_type,
                    sample,
                    analysis.case.key,
                )
                if track is None:
                    continue
                overlap = me_index._overlap_indices(
                    track,
                    analysis.case.key.position,
                )
                overlap = overlap[
                    track.support[overlap]
                    >= args.min_window_support
                ]
                if len(overlap) == 0:
                    continue
                values.append(
                    weighted_quantile(
                        track.entropy[overlap],
                        track.support[overlap],
                        0.50,
                    )
                )
            per_control_rows.append(
                {
                    "mod_type": mod_type,
                    "control_label": sample,
                    "case_sites_with_score": len(values),
                    "median_case_site_ME": (
                        float(np.median(values))
                        if values
                        else float("nan")
                    ),
                    "formal_replicate_inference": (
                        "not_performed; supply independent controls and use "
                        "a replicate-level model"
                    ),
                }
            )
    pd.DataFrame(per_control_rows).to_csv(
        output_dir / "results" / "control_track_summaries.tsv",
        sep="\t",
        index=False,
    )

    # Score table for reproducibility.
    score_rows = []
    for (mod_type, key, threshold), score in score_cache.items():
        if score is None:
            continue
        score_rows.append(
            {
                "mod_type": mod_type,
                "chrom": key.chrom,
                "start": key.start,
                "end": key.end,
                "strand": key.strand,
                "support_threshold": threshold,
                "strand_mode": score.strand_mode,
                "site_me": score.site_me,
                "site_q75": score.site_q75,
                "site_support": score.site_support,
                "overlap_windows": score.overlap_windows,
                "domain_250": score.domain_250,
                "domain_500": score.domain_500,
                "domain_1000": score.domain_1000,
                "context_500_excluding_site": (
                    score.context_500_excluding_site
                ),
                "local_prominence_500": score.local_prominence_500,
                "local_percentile_500": score.local_percentile_500,
                "local_context_windows_500": (
                    score.local_context_windows_500
                ),
                "local_window_median": score.local_window_median,
                "local_window_prominence": (
                    score.local_window_prominence
                ),
            }
        )
    pd.DataFrame(score_rows).to_csv(
        output_dir / "data" / "all_scored_sites.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    plot_case_percentiles(
        case_level,
        output_dir / "figures" / "placeholder.png",
    )
    write_summary(
        output_dir,
        primary,
        case_qc,
        me_qc,
        metadata_qc,
        me_inputs,
        provenance_qc,
        direction_frame,
        treatment_frame,
        sensitivity_frame,
        continuous_results,
    )

    # Hypothesis-first report requested for manuscript interpretation.
    write_hypothesis_results_report(
        output_dir,
        strict_matched=primary,
        continuous=continuous_results,
        binary=binary_results,
        me_differences=case_me_differences,
        direct=direct_peak_contrasts,
        union_sites=union_site_data,
    )

    hypothesis_registry = pd.DataFrame(
        [
            {
                "hypothesis_id": "H1",
                "plain_language_question": (
                    "Does higher pre-treatment site-overlapping ME predict "
                    "a larger later absolute PL change?"
                ),
                "primary_result_file": "continuous_response_models.tsv",
                "endpoint": "site_overlapping_ME",
            },
            {
                "hypothesis_id": "H2",
                "plain_language_question": (
                    "Does higher pre-treatment +/-500-bp neighbourhood ME "
                    "predict a larger later absolute PL change?"
                ),
                "primary_result_file": "continuous_response_models.tsv",
                "endpoint": "neighbourhood_500_ME",
            },
            {
                "hypothesis_id": "H3",
                "plain_language_question": (
                    "Does baseline site-overlapping ME predict whether a "
                    "site later crosses the significant response threshold?"
                ),
                "primary_result_file": "binary_significance_models.tsv",
                "endpoint": "site_overlapping_ME_per_0.1",
            },
            {
                "hypothesis_id": "H4",
                "plain_language_question": (
                    "Does baseline +/-500-bp neighbourhood ME predict "
                    "whether a site later crosses the response threshold?"
                ),
                "primary_result_file": "binary_significance_models.tsv",
                "endpoint": "neighbourhood_500_ME_per_0.1",
            },
            {
                "hypothesis_id": "H5",
                "plain_language_question": (
                    "How much higher is site-overlapping ME at later "
                    "threshold-positive sites than at the callable background?"
                ),
                "primary_result_file": "case_ME_difference_models.tsv",
                "endpoint": "site_overlapping_ME",
            },
            {
                "hypothesis_id": "H6",
                "plain_language_question": (
                    "How much higher is +/-500-bp neighbourhood ME at later "
                    "threshold-positive sites than at the callable background?"
                ),
                "primary_result_file": "case_ME_difference_models.tsv",
                "endpoint": "neighbourhood_500_ME",
            },
            {
                "hypothesis_id": "H7",
                "plain_language_question": (
                    "Are later threshold-positive sites focal ME peaks "
                    "relative to their own +/-500-bp context?"
                ),
                "primary_result_file": "direct_peak_contrasts.tsv",
                "endpoint": (
                    "site_overlapping_ME_vs_own_500bp_context"
                ),
            },
            {
                "hypothesis_id": "H8",
                "plain_language_question": (
                    "Are later threshold-positive sites high-ME positions "
                    "relative to the callable genome-wide background?"
                ),
                "primary_result_file": "direct_peak_contrasts.tsv",
                "endpoint": (
                    "site_overlapping_ME_vs_callable_genome_background"
                ),
            },
            {
                "hypothesis_id": "H9",
                "plain_language_question": (
                    "Do threshold-positive sites remain ME-enriched under "
                    "strict case-specific PL/depth/support matching?"
                ),
                "primary_result_file": "primary_tests.tsv",
                "endpoint": "strict_matched_local_and_global",
            },
        ]
    )
    hypothesis_registry.to_csv(
        output_dir / "results" / "hypothesis_registry.tsv",
        sep="\t",
        index=False,
    )

    write_methods_report(
        output_dir,
        args,
        me_qc,
    )

    # Environment versions.
    versions = pd.DataFrame(
        [
            {"component": "script", "version": VERSION},
            {"component": "python", "version": sys.version.replace("\n", " ")},
            {"component": "numpy", "version": np.__version__},
            {"component": "pandas", "version": pd.__version__},
            {"component": "matplotlib", "version": matplotlib.__version__},
        ]
    )
    versions.to_csv(
        output_dir / "software_versions.tsv",
        sep="\t",
        index=False,
    )

    connection.close()

    if not args.keep_work:
        for path in (
            database_path,
            Path(str(database_path) + "-wal"),
            Path(str(database_path) + "-shm"),
        ):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        try:
            work_dir = output_dir / "work"
            if work_dir.exists() and not any(work_dir.iterdir()):
                work_dir.rmdir()
        except OSError:
            pass

    zip_path = zip_results(output_dir)
    LOG.info("Analysis complete: %s", output_dir)
    LOG.info("Upload-ready archive: %s", zip_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
