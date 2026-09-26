#!/usr/bin/env python3
import argparse, collections, csv, sys, os
from statistics import median
from typing import Dict, Tuple, Set, List, Optional
import pysam

# ----------------------------- CLI -----------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        "Call Tn insertion sites from overhang→genome BAM, "
        "then compute local occupancy using an ALL-reads genome BAM, "
        "and emit extensive verification diagnostics to stderr."
    )
    ap.add_argument("--bam_over", required=True,
                    help="BAM: overhangs mapped to Tn-free genome (from overhang extraction step)")
    ap.add_argument("--bam_all", required=True,
                    help="BAM: ALL reads mapped to ALL contigs (Tn-free genome recommended; include Tn contig for global Tn-free estimation)")
    ap.add_argument("--meta", required=True,
                    help="Metadata TSV from overhang extraction (must include columns: oh_id, read, tn_strand)")
    ap.add_argument("--out", required=True, help="Output TSV")

    ap.add_argument("--min-mapq-genome", type=int, default=20,
                    help="Minimum MAPQ for overhang→genome alignments")
    ap.add_argument("--min-anchor-bp", type=int, default=60,
                    help="Minimum contiguous match near overhang 5' end (genome side)")
    ap.add_argument("--cluster", type=int, default=2,
                    help="Cluster window (bp) when collapsing sites")
    ap.add_argument("--min-support", type=int, default=3,
                    help="Filter sites with < this many supporting reads")

    # Circular merge implied by presence of --fai (no explicit --circular flag)
    ap.add_argument("--fai", help="Reference .fai (if provided, enables circular-origin merge using lengths)")

    # Local occupancy window for denominator from --bam_all
    ap.add_argument("--occupancy-window", type=int, default=200,
                    help="+/- bp window around breakpoint for local coverage counting (default 200)")

    # Diagnostics verbosity
    ap.add_argument("--diag-top-n", type=int, default=10,
                    help="Print detailed diagnostics for the top-N sites by support (default 10)")

    # NEW: global Tn-free prevalence (from combined-reference BAM)
    ap.add_argument("--tn-contig", help="Exact contig name of the Tn in --bam_all header to enable global Tn-free estimation")
    ap.add_argument("--tn-min-mapq", type=int, default=20,
                    help="Minimum MAPQ to count a Tn hit in global prevalence (default 20)")
    ap.add_argument("--tn-min-align-bp", type=int, default=100,
                    help="Minimum reference-aligned bp on Tn contig to count a Tn hit (default 100)")
    ap.add_argument("--tn-summary-out", help="Optional path to write a 1-row TSV summary of global Tn/Tn-free prevalence")
    ap.add_argument("--tn-write-readlists-prefix",
                    help="If set, write <prefix>.tn_pos.txt and <prefix>.tn_neg.txt with read IDs (global prevalence)")

    return ap.parse_args()


# ----------------------------- Utils -----------------------------
def ref_aligned_len(cigartuples) -> int:
    """
    Reference-consuming length from CIGAR: M(0), D(2), N(3), =(7), X(8).
    """
    if not cigartuples:
        return 0
    return sum(l for op, l in cigartuples if op in (0, 2, 3, 7, 8))


def tn_free_global_from_combined_bam(
    bam: pysam.AlignmentFile,
    tn_contig: str,
    min_mapq: int = 20,
    min_tn_align_bp: int = 100,
    write_lists_prefix: Optional[str] = None
) -> Dict[str, float]:
    """
    Classify molecules (read names) as Tn+ if ANY segment (primary or supplementary)
    aligns to tn_contig with MAPQ>=min_mapq and ref-aligned length>=min_tn_align_bp.
    Denominator = number of unique PRIMARY reads observed (molecules).

    Returns dict with counts and percentages. Prints a concise audit to stderr.
    Optionally writes read ID lists if write_lists_prefix is provided.
    """
    if tn_contig not in bam.references:
        raise ValueError(f"Tn contig '{tn_contig}' not found in BAM header.")

    seen_primary: Set[str] = set()
    tn_pos: Set[str] = set()
    tn_hit_lengths: List[int] = []
    tn_hit_mapqs: List[int] = []
    total_records = 0

    for aln in bam.fetch(until_eof=True):
        total_records += 1
        if aln.is_unmapped:
            continue

        # track primary read names for denominator
        if not aln.is_secondary and not aln.is_supplementary:
            seen_primary.add(aln.query_name)

        # classify Tn+ if any (primary or supplementary) segment passes thresholds on tn_contig
        if aln.reference_name == tn_contig and aln.mapping_quality >= min_mapq:
            L = ref_aligned_len(aln.cigartuples)
            if L >= min_tn_align_bp:
                tn_pos.add(aln.query_name)
                tn_hit_lengths.append(L)
                tn_hit_mapqs.append(aln.mapping_quality)

    total_reads = len(seen_primary)
    tn_reads = len(tn_pos)
    tn_free_reads = max(0, total_reads - tn_reads)
    tn_pos_pct = (100.0 * tn_reads / total_reads) if total_reads > 0 else 0.0
    tn_free_pct = (100.0 * tn_free_reads / total_reads) if total_reads > 0 else 0.0

    # audit
    stderr("\n==== GLOBAL Tn-FREE PREVALENCE (combined-ref) ====")
    stderr(f"tn_contig              : {tn_contig}")
    stderr(f"primary reads (denom)  : {total_reads}")
    stderr(f"Tn+ reads (criteria)   : {tn_reads}")
    stderr(f"Tn-free reads          : {tn_free_reads}")
    stderr(f"Tn+ prevalence  (%)    : {tn_pos_pct:.6f}")
    stderr(f"Tn-free prevalence (%) : {tn_free_pct:.6f}")
    if tn_hit_lengths:
        Ls = sorted(tn_hit_lengths)
        MQs = sorted(tn_hit_mapqs)
        q = lambda vs, p: vs[int((len(vs)-1)*p)]
        stderr(f"Tn hit len (bp) q10/50/90 : {q(Ls,0.10)}/{q(Ls,0.50)}/{q(Ls,0.90)}")
        stderr(f"Tn hit MAPQ q10/50/90     : {q(MQs,0.10)}/{q(MQs,0.50)}/{q(MQs,0.90)}")
    else:
        stderr("No Tn hits passed thresholds; all molecules classified as Tn-free.")

    # optional read lists
    if write_lists_prefix:
        with open(f"{write_lists_prefix}.tn_pos.txt","w") as outp:
            for r in sorted(tn_pos): outp.write(r+"\n")
        with open(f"{write_lists_prefix}.tn_neg.txt","w") as outn:
            for r in sorted(seen_primary - tn_pos): outn.write(r+"\n")
        stderr(f"[INFO] Wrote read lists: {write_lists_prefix}.tn_pos.txt / .tn_neg.txt")

    return {
        "primary_reads": float(total_reads),
        "tn_pos_reads": float(tn_reads),
        "tn_free_reads": float(tn_free_reads),
        "tn_pos_pct": tn_pos_pct,
        "tn_free_pct": tn_free_pct,
    }

def stderr(s: str):
    print(s, file=sys.stderr)

def contiguous_match_near_end_ct(ct, near_left: bool) -> int:
    """
    Return contiguous matched length (M/= /X) adjacent to the read end that is
    proximal to the junction. Accepts pysam cigartuples; counts S/H as soft/hard clips.
    """
    if not ct:
        return 0
    ops = list(ct)
    if near_left:
        i = 0
        if i < len(ops) and ops[i][0] in (4, 5):  # S or H
            i += 1
        return ops[i][1] if i < len(ops) and ops[i][0] in (0, 7, 8) else 0  # M,=,X
    else:
        i = len(ops) - 1
        if i >= 0 and ops[i][0] in (4, 5):
            i -= 1
        return ops[i][1] if i >= 0 and ops[i][0] in (0, 7, 8) else 0

def load_meta(path: str) -> Dict[str, Dict[str, str]]:
    """
    Load metadata keyed by overhang id (oh_id). Must have at least:
      oh_id, read, tn_strand
    """
    meta: Dict[str, Dict[str, str]] = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            meta[row["oh_id"]] = row
    return meta

def has_bai(bam: pysam.AlignmentFile) -> bool:
    try:
        return bam.has_index()
    except Exception:
        try:
            bam.check_index()
            return True
        except Exception:
            return False

def local_unique_reads(bam: pysam.AlignmentFile, chrom: str, start0: int, end0: int) -> int:
    """
    Count unique primary, mapped read names overlapping [start0, end0).
    Uses indexed fetch when available; otherwise falls back to a linear scan.
    """
    names: Set[str] = set()
    if has_bai(bam):
        for a in bam.fetch(chrom, max(0, start0), end0):
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            names.add(a.query_name)
        return len(names)
    # Fallback: linear scan (slower, but robust when no index present)
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        if a.reference_name != chrom:
            continue
        if not (a.reference_end <= start0 or a.reference_start >= end0):
            names.add(a.query_name)
    return len(names)

def quantiles(values: List[int], qs=(0.1, 0.25, 0.5, 0.75, 0.9)) -> Dict[float, float]:
    if not values:
        return {q: 0.0 for q in qs}
    vs = sorted(values)
    out={}
    n=len(vs)
    for q in qs:
        idx = q*(n-1)
        lo = int(idx); hi = min(n-1, lo+1)
        frac = idx - lo
        out[q] = vs[lo]*(1-frac) + vs[hi]*frac
    return out

def bam_quick_stats(bam: pysam.AlignmentFile, min_mapq: int) -> Dict[str, int]:
    total=0; mapped=0; primary=0; primary_mq=0
    mq_list=[]
    for a in bam.fetch(until_eof=True):
        total += 1
        if not a.is_unmapped:
            mapped += 1
            if not a.is_secondary and not a.is_supplementary:
                primary += 1
                if a.mapping_quality >= min_mapq:
                    primary_mq += 1
                    mq_list.append(a.mapping_quality)
    return {
        "total": total,
        "mapped": mapped,
        "primary": primary,
        "primary_mq_ge": primary_mq,
        "mq_min": min(mq_list) if mq_list else -1,
        "mq_max": max(mq_list) if mq_list else -1,
        "mq_median": int(median(mq_list)) if mq_list else -1
    }

# ----------------------------- Main -----------------------------
def main():
    a = parse_args()

    # ---- Parameter echo
    stderr("==== PARAMETERS ====")
    for k,v in [
        ("bam_over", a.bam_over),
        ("bam_all", a.bam_all),
        ("meta", a.meta),
        ("out", a.out),
        ("min_mapq_genome", a.min_mapq_genome),
        ("min_anchor_bp", a.min_anchor_bp),
        ("cluster", a.cluster),
        ("min_support", a.min_support),
        ("fai", a.fai if a.fai else "None"),
        ("occupancy_window", a.occupancy_window),
        ("diag_top_n", a.diag_top_n),
    ]:
        stderr(f"{k:20s}: {v}")

    # ---- File presence
    for p in [a.bam_over, a.bam_all, a.meta]:
        if not os.path.exists(p):
            sys.exit(f"[ERROR] Missing file: {p}")

    # Open inputs
    meta = load_meta(a.meta)
    if not meta:
        sys.exit("[ERROR] Meta is empty or missing required columns (oh_id, read, tn_strand).")

    bam_over = pysam.AlignmentFile(a.bam_over, "rb")
    bam_all  = pysam.AlignmentFile(a.bam_all,  "rb")

    # ---- Index status
    stderr("\n==== INDEX STATUS ====")
    stderr(f"bam_over indexed   : {has_bai(bam_over)}")
    stderr(f"bam_all  indexed   : {has_bai(bam_all)}")

    # ---- Header / contig consistency
    over_refs = list(bam_over.references)
    all_refs  = list(bam_all.references)
    shared    = sorted(set(over_refs) & set(all_refs))
    only_over = sorted(set(over_refs) - set(all_refs))
    only_all  = sorted(set(all_refs)  - set(over_refs))
    stderr("\n==== CONTIG CONSISTENCY ====")
    stderr(f"overhang BAM contigs: {len(over_refs)}")
    stderr(f"all-reads BAM contigs: {len(all_refs)}")
    stderr(f"shared contigs      : {len(shared)}")
    if only_over:
        stderr(f"[WARN] contigs only in bam_over (first 5): {only_over[:5]}")
    if only_all:
        stderr(f"[WARN] contigs only in bam_all (first 5): {only_all[:5]}")

    # Reference lengths for circular merge (enabled iff --fai provided)
    ref_lens: Dict[str, int] = {r: L for r, L in zip(bam_over.references, bam_over.lengths)}
    circular = False
    if a.fai:
        circular = True
        with open(a.fai) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    ref_lens[parts[0]] = int(parts[1])
        # Optional: check fai contigs
        fai_names = []
        with open(a.fai) as fh:
            for line in fh:
                fai_names.append(line.split("\t",1)[0])
        miss_in_over = sorted(set(fai_names) - set(over_refs))
        if miss_in_over:
            stderr(f"[WARN] contigs in FAI but NOT in bam_over header (first 5): {miss_in_over[:5]}")

    # ---- Quick BAM stats
    stderr("\n==== BAM SUMMARIES ====")
    s_over = bam_quick_stats(pysam.AlignmentFile(a.bam_over, "rb"), a.min_mapq_genome)
    s_all  = bam_quick_stats(pysam.AlignmentFile(a.bam_all,  "rb"), a.min_mapq_genome)
    stderr(f"bam_over total     : {s_over['total']}")
    stderr(f"bam_over mapped    : {s_over['mapped']} (primary: {s_over['primary']}, primary MAPQ≥{a.min_mapq_genome}: {s_over['primary_mq_ge']})")
    stderr(f"bam_over MAPQ      : min={s_over['mq_min']} median={s_over['mq_median']} max={s_over['mq_max']}")
    stderr(f"bam_all  total     : {s_all['total']}")
    stderr(f"bam_all  mapped    : {s_all['mapped']} (primary: {s_all['primary']}, primary MAPQ≥{a.min_mapq_genome}: {s_all['primary_mq_ge']})")
    stderr(f"bam_all  MAPQ      : min={s_all['mq_min']} median={s_all['mq_median']} max={s_all['mq_max']}")


    # -------------------- Global Tn-free prevalence (optional) --------------------
    if a.tn_contig:
        try:
            # Use bam_all as the combined-reference BAM; it must include the Tn contig
            if a.tn_contig not in bam_all.references:
                stderr(f"\n[WARN] tn_contig '{a.tn_contig}' not present in --bam_all header; "
                       f"global Tn-free estimation skipped.")
            else:
                tn_stats = tn_free_global_from_combined_bam(
                    bam_all,
                    tn_contig=a.tn_contig,
                    min_mapq=a.tn_min_mapq,
                    min_tn_align_bp=a.tn_min_align_bp,
                    write_lists_prefix=a.tn_write_readlists_prefix
                )
                # Optional 1-row TSV summary
                if a.tn_summary_out:
                    with open(a.tn_summary_out, "w") as fh:
                        fh.write("\t".join([
                            "bam","tn_contig","min_mapq","min_tn_align_bp",
                            "primary_reads","tn_pos_reads","tn_free_reads","tn_pos_pct","tn_free_pct"
                        ]) + "\n")
                        fh.write("\t".join([
                            a.bam_all, a.tn_contig, str(a.tn_min_mapq), str(a.tn_min_align_bp),
                            str(int(tn_stats["primary_reads"])),
                            str(int(tn_stats["tn_pos_reads"])),
                            str(int(tn_stats["tn_free_reads"])),
                            f"{tn_stats['tn_pos_pct']:.6f}",
                            f"{tn_stats['tn_free_pct']:.6f}",
                        ]) + "\n")
                        stderr(f"[INFO] Wrote global Tn-free summary TSV: {a.tn_summary_out}")
        except Exception as e:
            stderr(f"[ERROR] Global Tn-free estimation failed: {e}")
    else:
        stderr("\n==== GLOBAL Tn-FREE PREVALENCE ====")
        stderr("tn_contig not provided; skipping global Tn-free estimation.")

    # -------------------- Discover sites from overhang→genome --------------------
    # Keyed by (chrom, pos0, genome_strand, tn_strand)
    sites: Dict[Tuple[str, int, str, str], Set[str]] = collections.defaultdict(set)
    site_anchor: Dict[Tuple[str, int, str, str], List[int]] = collections.defaultdict(list)
    anchor_lengths_all: List[int] = []

    for aln in pysam.AlignmentFile(a.bam_over, "rb").fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary:
            continue
        if aln.mapping_quality < a.min_mapq_genome:
            continue

        qname = aln.query_name  # equals oh_id from overhang FASTA header
        info = meta.get(qname)
        if info is None:
            continue

        tn_strand = info["tn_strand"]  # '+' or '-'
        genome_strand = "-" if aln.is_reverse else "+"

        # Choose junction-proximal reference coordinate and anchor length
        if genome_strand == "+":
            pos0 = aln.reference_start
            anchor_bp = contiguous_match_near_end_ct(aln.cigartuples, near_left=True)
        else:
            pos0 = aln.reference_end
            anchor_bp = contiguous_match_near_end_ct(aln.cigartuples, near_left=False)

        if anchor_bp < a.min_anchor_bp:
            continue

        anchor_lengths_all.append(anchor_bp)
        key = (aln.reference_name, pos0, genome_strand, tn_strand)
        sites[key].add(info["read"])
        site_anchor[key].append(anchor_bp)

    stderr("\n==== OVERHANG DISCOVERY ====")
    stderr(f"raw junction keys       : {len(sites)}")
    stderr(f"anchor length (bp) q10/q25/median/q75/q90: " +
           "/".join(f"{quantiles(anchor_lengths_all)[q]:.1f}" for q in (0.1,0.25,0.5,0.75,0.9)) if anchor_lengths_all else "NA")

    # -------------------- Cluster within window --------------------
    buckets: Dict[Tuple[str, str, str], List[Tuple[int, Set[str], List[int]]]] = collections.defaultdict(list)
    for (chrom, pos0, gs, ts), rset in sites.items():
        buckets[(chrom, gs, ts)].append((pos0, rset, site_anchor[(chrom, pos0, gs, ts)]))

    clustered: List[Tuple[str, int, str, str, Set[str], List[int]]] = []
    for (chrom, gs, ts), arr in buckets.items():
        arr.sort(key=lambda x: x[0])
        cur_pos: Optional[int] = None
        cur_reads: Set[str] = set()
        cur_anchors: List[int] = []
        for pos0, rset, anchors in arr:
            if cur_pos is None or abs(pos0 - cur_pos) > a.cluster:
                if cur_pos is not None:
                    clustered.append((chrom, cur_pos, gs, ts, cur_reads, cur_anchors))
                cur_pos = pos0
                cur_reads = set(rset)
                cur_anchors = list(anchors)
            else:
                cur_reads |= rset
                cur_anchors += anchors
                cur_pos = (cur_pos + pos0) // 2
        if cur_pos is not None:
            clustered.append((chrom, cur_pos, gs, ts, cur_reads, cur_anchors))

    stderr("\n==== CLUSTERING ====")
    stderr(f"cluster window (bp)      : {a.cluster}")
    stderr(f"clustered sites (pre-circ): {len(clustered)}")

    # -------------------- Optional circular merge across origin --------------------
    if a.fai:
        bychr: Dict[str, List[Tuple[str, int, str, str, Set[str], List[int]]]] = collections.defaultdict(list)
        for rec in clustered:
            bychr[rec[0]].append(rec)
        merged: List[Tuple[str, int, str, str, Set[str], List[int]]] = []
        merges_done = 0
        for chrom, items in bychr.items():
            L = ref_lens.get(chrom, 0)
            if not L or len(items) < 2:
                merged += items
                continue
            items.sort(key=lambda x: x[1])
            first, last = items[0], items[-1]
            if first[1] <= a.cluster and (L - last[1]) <= a.cluster and (first[2], first[3]) == (last[2], last[3]):
                merged_reads   = set(first[4]) | set(last[4])
                merged_anchors = first[5] + last[5]
                merged.append((chrom, 0, first[2], first[3], merged_reads, merged_anchors))
                items = items[1:-1]
                merges_done += 1
            merged += items
        clustered = merged
        stderr("\n==== CIRCULAR MERGE ====")
        stderr(f"FAI provided; merges across origin performed: {merges_done}")
    else:
        stderr("\n==== CIRCULAR MERGE ====")
        stderr("FAI not provided; circular merge skipped.")

    # -------------------- Filter by minimum support --------------------
    pre_filter = len(clustered)
    clustered = [rec for rec in clustered if len(rec[4]) >= a.min_support]
    stderr("\n==== FILTERING ====")
    stderr(f"min_support             : {a.min_support}")
    stderr(f"sites kept              : {len(clustered)} (dropped {pre_filter - len(clustered)})")

    # -------------------- Denominator for RPM/prevalence (bin sum) --------------------
    denom_bin = float(sum(len(rec[4]) for rec in clustered) or 1.0)
    stderr("\n==== DENOMINATOR (GLOBAL BIN) ====")
    stderr(f"bin-sum denominator     : {denom_bin:.0f} (sum of supports across kept sites)")

    # -------------------- Output with local occupancy from ALL-reads BAM --------------------
    occ_win = int(a.occupancy_window)
    stderr("\n==== LOCAL OCCUPANCY ====")
    stderr(f"occupancy window ±bp    : {occ_win}")

    # Write and simultaneously collect top-N diagnostics
    rows_for_diag: List[Tuple[int, int, float, str, int, str, str]] = []  # (support, local_cov, occ_pct, chrom, pos1, gs, ts)

    with open(a.out, "w") as out:
        out.write("\t".join([
            "chrom", "pos_1based", "genome_strand", "tn_strand",
            "support_reads", "RPM", "prevalence_pct",
            "median_anchor_bp", "read_names",
            "local_coverage_reads", "occupancy_local_pct"
        ]) + "\n")

        for chrom, pos0, gs, ts, reads, anchors in clustered:
            support = len(reads)
            rpm = 1e6 * support / denom_bin
            prev = 100.0 * support / denom_bin
            med_anchor = int(median(anchors)) if anchors else 0
            preview = ",".join(list(reads)[:5]) + (";..." if len(reads) > 5 else "")

            # Local denominator from ALL-reads BAM (unique primary molecules overlapping a window)
            start = max(0, pos0 - occ_win)
            end   = pos0 + occ_win + 1
            local_cov = local_unique_reads(bam_all, chrom, start, end)
            occ_pct = (100.0 * support / local_cov) if local_cov > 0 else 0.0

            # Safety flags
            if local_cov == 0:
                stderr(f"[WARN ZERO COV] {chrom}:{pos0+1} support={support} local_cov=0 (check contigs & reference)")
            if local_cov > 0 and local_cov < support:
                stderr(f"[WARN COV<NUM] {chrom}:{pos0+1} support={support} local_cov={local_cov} (window too small?)")

            out.write(
                f"{chrom}\t{pos0+1}\t{gs}\t{ts}\t{support}\t{rpm:.2f}\t{prev:.3f}\t"
                f"{med_anchor}\t{preview}\t{local_cov}\t{occ_pct:.6f}\n"
            )

            rows_for_diag.append((support, local_cov, occ_pct, chrom, pos0+1, gs, ts))

    # -------------------- Top-N site diagnostics --------------------
    rows_for_diag.sort(key=lambda x: x[0], reverse=True)
    topN = rows_for_diag[:max(0, a.diag_top_n)]
    stderr("\n==== TOP-SITE DIAGNOSTICS ====")
    if not topN:
        stderr("No sites kept after filtering.")
    else:
        stderr("support\tlocal_cov\tocc_pct\tchrom\tpos_1based\tg_strand\ttn_strand")
        for (support, local_cov, occ_pct, chrom, pos1, gs, ts) in topN:
            stderr(f"{support}\t{local_cov}\t{occ_pct:.6f}\t{chrom}\t{pos1}\t{gs}\t{ts}")

    # -------------------- Final summary --------------------
    stderr("\n==== SUMMARY ====")
    stderr(f"sites written           : {len(rows_for_diag)}")
    stderr(f"output                  : {a.out}")
    stderr("Done.")

if __name__ == "__main__":
    main()
