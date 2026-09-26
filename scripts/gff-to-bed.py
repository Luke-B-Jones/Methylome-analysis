#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert GFF3 to BED6 using locus_tag as name (ignores ANONYMOUS)")
    parser.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()

    with open(args.input) as inf, open(args.output, "w") as outf:
        for line in inf:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = parts

            # convert to 0-based BED
            start0 = str(int(start) - 1)
            end1   = end

            # extract locus_tag
            locus = None
            for field in attrs.split(";"):
                if field.startswith("locus_tag="):
                    locus = field.split("=", 1)[1]
                    break

            if not locus:
                continue  # skip if no locus_tag

            # ignore ANONYMOUS genes
            if locus == "ANONYMOUS":
                continue

            # BED6: chrom, start, end, locus_tag, score=0, strand
            outf.write("\t".join([chrom, start0, end1, locus, "0", strand]) + "\n")

if __name__ == "__main__":
    main()
