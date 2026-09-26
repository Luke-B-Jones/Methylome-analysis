#!/usr/bin/env python3
import argparse, collections
from typing import List, Tuple, Optional
import pysam

def parse_args():
    ap = argparse.ArgumentParser("Extract soft-clipped overhangs from reads mapped to Tn")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--tn-contig", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-meta", required=True)
    ap.add_argument("--min-mapq-tn", type=int, default=20)
    ap.add_argument("--min-overhang", type=int, default=40)
    ap.add_argument("--tn-end-window", type=int, default=40)
    return ap.parse_args()

def clip_lengths(ct: Optional[List[Tuple[int,int]]]) -> Tuple[int,int]:
    if not ct: return (0,0)
    L = ct[0][1] if ct[0][0] in (4,5) else 0
    R = ct[-1][1] if ct[-1][0] in (4,5) else 0
    return L,R

def ref_span(ct: Optional[List[Tuple[int,int]]]) -> int:
    if not ct: return 0
    s=0
    for op,len_ in ct:
        if op in (0,2,3,7,8): # M,D,N,=,X consume ref
            s+=len_
    return s

def main():
    a = parse_args()
    bam = pysam.AlignmentFile(a.bam, "rb")
    tn_len = dict(zip(bam.references, bam.lengths)).get(a.tn_contig, 0)
    if tn_len == 0:
        raise SystemExit(f"Tn contig {a.tn_contig} not found in BAM header")

    # choose best Tn alignment per (read, side) to avoid duplicates
    best = {}  # key: (qname, side) -> (aln, overhang_len)
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary: continue
        if aln.reference_name != a.tn_contig: continue
        if aln.mapping_quality < a.min_mapq_tn: continue
        if aln.cigartuples is None or aln.query_sequence is None: continue
        L,R = clip_lengths(aln.cigartuples)
        if L < a.min_overhang and R < a.min_overhang:  # need at least one usable overhang
            continue
        if L >= a.min_overhang:
            key = (aln.query_name, "L")
            cand = (aln, L)
            if key not in best or L > best[key][1]:
                best[key] = cand
        if R >= a.min_overhang:
            key = (aln.query_name, "R")
            cand = (aln, R)
            if key not in best or R > best[key][1]:
                best[key] = cand

    # write fasta + meta
    n=0
    with open(a.out_fasta, "w") as fa, open(a.out_meta, "w") as mt:
        mt.write("oh_id\tread\tside\ttn_strand\ttn_start0\ttn_end0\ttn_end_dist\toverhang_len\n")
        for (qname, side),(aln, olen) in best.items():
            seq = aln.query_sequence
            if seq is None: continue
            # extract overhang bases (soft-clipped)
            ct = aln.cigartuples
            L,R = clip_lengths(ct)
            if side == "L":
                oh = seq[:L]
            else:
                oh = seq[-R:]
            if len(oh) < a.min_overhang: continue

            tn_start0, tn_end0 = aln.reference_start, aln.reference_end
            # distance to nearest Tn TERMINUS (0 or tn_len)
            dist = min(tn_start0, tn_len - tn_end0)

            # require overhang to be from near a Tn end (optional but strong)
            if dist > a.tn_end_window:
                continue

            tn_strand = "-" if aln.is_reverse else "+"
            oh_id = f"{qname}|{side}"
            fa.write(f">{oh_id}\n{oh}\n")
            mt.write(f"{oh_id}\t{qname}\t{side}\t{tn_strand}\t{tn_start0}\t{tn_end0}\t{dist}\t{len(oh)}\n")
            n+=1
    print(f"[extract_overhangs] wrote {n} overhangs")

if __name__ == "__main__":
    main()
