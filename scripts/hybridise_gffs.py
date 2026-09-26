#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
import copy

def parse_attributes(attr_str):
    """Convert a semicolon‐separated attributes string into a dictionary."""
    attrs = {}
    if attr_str.strip() == "":
        return attrs
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
        else:
            attrs[part.strip()] = ""
    return attrs

def format_attributes(attrs):
    """Convert an attributes dictionary back to a semicolon‐separated string."""
    return ";".join(f"{k}={v}" for k, v in attrs.items())

def parse_gff(file_path):
    """
    Parse a GFF file into a list of feature dictionaries.
    Supports standard 9-column files as well as minimal 6-column files.
    For 6-column files, defaults are:
      source: "RNA" (overridden below),
      phase: ".",
      attributes: ""
    """
    features = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 9:
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    sys.stderr.write(f"Warning: Invalid numeric values in line: {line}\n")
                    continue
                feature = {
                    'seqid': parts[0],
                    'source': parts[1],
                    'type': parts[2],
                    'start': start,
                    'end': end,
                    'score': parts[5],
                    'strand': parts[6],
                    'phase': parts[7],
                    'attributes': parts[8]
                }
                features.append(feature)
            elif len(parts) == 6:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    sys.stderr.write(f"Warning: Invalid numeric values in line: {line}\n")
                    continue
                feature = {
                    'seqid': parts[0],
                    'source': "RNA",
                    'type': parts[3],
                    'start': start,
                    'end': end,
                    'score': parts[4],
                    'strand': parts[5],
                    'phase': ".",
                    'attributes': ""
                }
                features.append(feature)
            else:
                sys.stderr.write(f"Warning: Skipping malformed line (unexpected number of columns): {line}\n")
    return features

def write_gff(features, out_file):
    """Write a list of feature dictionaries to a GFF file."""
    with open(out_file, 'w') as f:
        f.write("##gff-version 3\n")
        for feat in features:
            line = "\t".join([
                feat['seqid'],
                feat['source'],
                feat['type'],
                str(feat['start']),
                str(feat['end']),
                feat['score'],
                feat['strand'],
                feat['phase'],
                feat['attributes']
            ])
            f.write(line + "\n")

def deduplicate_features(features):
    """
    Deduplicate features based on (seqid, start, end, strand).
    If duplicates are found, the first encountered is retained.
    """
    dedup = {}
    for feat in features:
        key = (feat['seqid'], feat['start'], feat['end'], feat['strand'])
        if key in dedup:
            continue
        dedup[key] = feat
    return list(dedup.values())

def update_feature_names(features, file_type):
    """
    Update each feature:
      - Override the source field with file_type ("PGAP" or "RNA").
      - For PGAP features: use the existing locus_tag to set the ID and warn if missing.
      - For RNA features: set the type field to '.' and move the original type into attributes
        (as Name, locus_tag, and ID) without warning.
    """
    for feat in features:
        feat['source'] = file_type
        if file_type == "RNA":
            original_type = feat['type']
            feat['type'] = "."
            attrs = parse_attributes(feat['attributes'])
            if 'Name' not in attrs:
                attrs['Name'] = original_type
            if 'locus_tag' not in attrs:
                attrs['locus_tag'] = original_type
            if 'ID' not in attrs:
                attrs['ID'] = original_type
            feat['attributes'] = format_attributes(attrs)
        else:  # PGAP features
            attrs = parse_attributes(feat['attributes'])
            if 'locus_tag' not in attrs:
                sys.stderr.write(f"Warning: locus_tag not found for PGAP feature at {feat['seqid']}:{feat['start']}-{feat['end']}\n")
            else:
                locus = attrs['locus_tag']
                attrs['ID'] = locus
            feat['attributes'] = format_attributes(attrs)
    return features

def intervals_overlap(s1, e1, s2, e2):
    """Return True if two intervals overlap."""
    return s1 <= e2 and s2 <= e1

def add_rna_features(rna_features, pgap_features):
    """
    For each RNA feature (from the RNA file), check if it overlaps any PGAP feature on the same seqid and strand.
    Only add the RNA feature if no overlap is found.
    Additionally, if the same RNA feature (same seqid, start, end, strand, original type) appears more than once,
    then add an extra copy with an attribute 'replicated' set to "*bacteria*".
    """
    grouped = defaultdict(list)
    for rna in rna_features:
        attrs = parse_attributes(rna['attributes'])
        name_val = attrs.get("Name", rna['type'])
        key = (rna['seqid'], rna['start'], rna['end'], rna['strand'], name_val)
        grouped[key].append(rna)
    
    added = 0
    for key, group in grouped.items():
        rep_rna = group[0]
        overlap_found = False
        for feat in pgap_features:
            # Check for same seqid and strand before testing coordinate overlap.
            if rep_rna['seqid'] != feat['seqid'] or rep_rna['strand'] != feat['strand']:
                continue
            if intervals_overlap(rep_rna['start'], rep_rna['end'], feat['start'], feat['end']):
                overlap_found = True
                break
        if not overlap_found:
            rep_rna['source'] = "RNA"
            attrs = parse_attributes(rep_rna['attributes'])
            if 'locus_tag' in attrs:
                locus = attrs['locus_tag']
            else:
                locus = key[4]
            attrs['ID'] = locus
            rep_rna['attributes'] = format_attributes(attrs)
            pgap_features.append(rep_rna)
            added += 1
            if len(group) > 1:
                rep2 = copy.deepcopy(rep_rna)
                attrs2 = parse_attributes(rep2['attributes'])
                attrs2['replicated'] = "*bacteria*"
                rep2['attributes'] = format_attributes(attrs2)
                pgap_features.append(rep2)
                added += 1
    sys.stderr.write(f"Internal check: Added {added} RNA feature(s) (including duplicate handling) not covered by PGAP.\n")
    return pgap_features

def main():
    parser = argparse.ArgumentParser(
        description="Merge and process PGAP and RNA GFF files: filter out unwanted homology lines, deduplicate by coordinates, "
                    "update attributes using locus_tag, override column 2 with PGAP or RNA, and for RNA features set type to '.' "
                    "with original type stored in Name, locus_tag, and ID."
    )
    parser.add_argument("--RNA", required=True, help="Path to the RNA GFF file (ncRNA annotations from cmsearch; may be 6-column or 9-column)")
    parser.add_argument("--PGAP", required=True, help="Path to the PGAP GFF file")
    parser.add_argument("-o", required=True, help="Output file path for the merged GFF")
    args = parser.parse_args()

    # Process PGAP file first.
    pgap_features = parse_gff(args.PGAP)
    sys.stderr.write(f"Internal check: Parsed {len(pgap_features)} features from PGAP file.\n")
    pgap_features = [
        feat for feat in pgap_features
        if feat['type'].lower() != "homology" and "protein homology" not in feat['source'].lower()
    ]
    sys.stderr.write(f"Internal check: {len(pgap_features)} PGAP features remain after filtering homology lines.\n")
    dedup_pgaps = deduplicate_features(pgap_features)
    sys.stderr.write(f"Internal check: {len(dedup_pgaps)} PGAP features remain after deduplication.\n")
    dedup_pgaps = update_feature_names(dedup_pgaps, "PGAP")

    # Process RNA file.
    rna_features = parse_gff(args.RNA)
    sys.stderr.write(f"Internal check: Parsed {len(rna_features)} features from RNA file.\n")
    rna_features = update_feature_names(rna_features, "RNA")

    # Merge RNA features that do not overlap PGAP features.
    merged_features = add_rna_features(rna_features, dedup_pgaps)
    sys.stderr.write(f"Internal check: Total merged features count is {len(merged_features)}.\n")
    
    # Sort final features by contig (seqid) and then by start coordinate.
    merged_features.sort(key=lambda x: (x['seqid'], x['start']))
    
    # Final pass: For RNA features only, update Name, ID, and locus_tag to be unique by appending a numeric suffix.
    rna_groups = defaultdict(list)
    for feat in merged_features:
        if feat['source'] == "RNA":
            attrs = parse_attributes(feat['attributes'])
            name = attrs.get("Name", "")
            rna_groups[name].append(feat)
    
    for name, feats in rna_groups.items():
        if len(feats) > 0:
            for i, feat in enumerate(feats, start=1):
                attrs = parse_attributes(feat['attributes'])
                new_val = f"{name}_{i}"
                attrs["Name"] = new_val
                attrs["ID"] = new_val
                attrs["locus_tag"] = new_val
                feat['attributes'] = format_attributes(attrs)
    
    # Write merged output.
    write_gff(merged_features, args.o)
    sys.stderr.write(f"Output written to {args.o}\n")

if __name__ == "__main__":
    main()
