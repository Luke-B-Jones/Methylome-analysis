import re
import argparse

def extract_motifs(blast_file, output_tsv):
    with open(blast_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        outfile.write("Query ID\tBinding Motif\n")

        query_motifs = {}

        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            query_id = parts[0].replace("gnl|extdb|", "")
            subject_id = parts[1]

            # Extract motif
            motif_match = re.search(r"::([A-Z]+[A-Z0-9]*)::", subject_id)
            motif = motif_match.group(1) if motif_match else "None"

            if query_id in query_motifs:
                if motif != "None" and motif not in query_motifs[query_id]:
                    query_motifs[query_id].append(motif)
            else:
                query_motifs[query_id] = [motif] if motif != "None" else []

        for query_id, motifs in query_motifs.items():
            motif_str = ",".join(motifs) if motifs else "None"
            outfile.write(f"{query_id}\t{motif_str}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract binding motifs from BLAST output.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input BLAST output file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    args = parser.parse_args()

    extract_motifs(args.input, args.output)

if __name__ == "__main__":
    main()
