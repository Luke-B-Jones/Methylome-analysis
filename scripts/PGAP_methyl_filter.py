import argparse
import re

def extract_locus_tags_and_functions(input_file, output_file):
    # Keywords that should match only at word boundaries
    boundary_keywords = [
        "restriction", "endonuclease", "modification",
        "type I", "type II", "subunit M", "subunit S", "subunit R",
        "Hsd", "RM system"
    ]
    # Keywords that can match inside compound words
    relaxed_keywords = [
        "methyltransfer", "methylase", "DNA methyl"
    ]

    boundary = r"(?:^|[\s\-])"

    # Compile combined regex
    pattern = re.compile(
        "|".join(
            [fr"{boundary}{kw}(?=$|[\s\-])" for kw in boundary_keywords] +
            [kw for kw in relaxed_keywords]
        ),
        re.IGNORECASE
    )

    locus_tag_pattern = re.compile(r"pgaptmp_\d+")
    results = []

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith(">") and pattern.search(line):
                tag_match = locus_tag_pattern.search(line)
                if tag_match:
                    locus_tag = tag_match.group()
                    function_desc = line[1:].strip()
                    results.append((locus_tag, function_desc))

    with open(output_file, 'w') as out:
        for tag, func in results:
            out.write(f"{tag}\t{func}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract locus tags and functional annotations for methylation-related proteins.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA (.faa) file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    args = parser.parse_args()
    extract_locus_tags_and_functions(args.input, args.output)
