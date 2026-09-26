#!/usr/bin/env python3
"""
Merge three TSV annotation files by locus tag into a single TSV table.

-r / --rebase    : rebase.tsv, with header "Query ID<TAB>Binding Motif"
- i / --interpro : interpro.tsv, with locus in col1 (maybe "gnl|extdb|<locus>") and description in col5
- p / --pgap      : pgap.tsv, with locus in col1 and description in col2
- o / --output    : out.tsv to write:
    locus_tag (col1), rebase binding motif (col2), all interpro descs joined by ';' (col3), pgap description (col4)
"""
import csv
import argparse
import sys

def parse_rebase(path):
    d = {}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        # skip header
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            locus, motif = row[0].strip(), row[1].strip()
            if locus:
                d[locus] = motif
    return d

def parse_interpro(path):
    d = {}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 5:
                continue
            raw = row[0].strip()
            # drop any gnl|extdb| prefix
            if '|' in raw:
                locus = raw.split('|')[-1]
            else:
                locus = raw
            desc = row[4].strip()
            if locus:
                d.setdefault(locus, []).append(desc)
    return d

def parse_pgap(path):
    d = {}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue
            locus, desc = row[0].strip(), row[1].strip()
            if locus:
                d[locus] = desc
    return d

def main():
    p = argparse.ArgumentParser(
        description="Merge rebase, interpro & pgap TSVs by locus tag")
    p.add_argument('-r', '--rebase',  required=True, help='rebase.tsv')
    p.add_argument('-i', '--interpro', required=True, help='interpro.tsv')
    p.add_argument('-p', '--pgap',    required=True, help='pgap.tsv')
    p.add_argument('-o', '--output',  required=True, help='out.tsv')
    args = p.parse_args()

    rebase    = parse_rebase(args.rebase)
    interpro  = parse_interpro(args.interpro)
    pgap      = parse_pgap(args.pgap)

    all_loci = set(rebase) | set(interpro) | set(pgap)
    if not all_loci:
        sys.exit("No locus tags found in any input files.")

    with open(args.output, 'w', newline='', encoding='utf-8') as fout:
        writer = csv.writer(fout, delimiter='\t')
        # header
        writer.writerow([
            "LocusTag",
            "Rebase_BindingMotif",
            "InterPro_Descriptions",
            "PGAP_Description"
        ])
        for locus in sorted(all_loci):
            motif    = rebase.get(locus, "")
            ip_descs = interpro.get(locus, [])
            ip_field = ";".join(ip_descs) if ip_descs else ""
            pgap_desc= pgap.get(locus, "")
            writer.writerow([locus, motif, ip_field, pgap_desc])

    print(f"Done: merged {len(all_loci)} loci into {args.output}")

if __name__ == '__main__':
    main()
