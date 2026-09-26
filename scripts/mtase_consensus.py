#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def get_search_terms(args_type):
    """
    Returns one of two search-term lists, based on the requested --type:
    - MTASE: Methyltransferase-only and RM system components
    - ALL:   Methyltransferase + DNA-binding/regulator terms
    """

    # 1) Exclusively Methyltransferases and RM components (all types)
    search_terms_mtase = [
        # General MTase terms
        "methyltransferase", "methylase", "DNA methyltransferase", "RNA methyltransferase",
        "SAM-dependent", "S-adenosylmethionine", "N6-adenine", "C5-cytosine", "DNA adenine methylase", "Dam", "Dcm",

        # RM system types
        "restriction-modification", "type I restriction", "type II restriction",
        "type III restriction", "type IV restriction",

        # RM subunits and domain types
        "restriction endonuclease", "endonuclease subunit", "modification methyltransferase",
        "specificity subunit", "HsdR", "HsdM", "HsdS",

        # Domain/family identifiers
        "IPR003356",  # N6-adenine-specific DNA MTase
        "IPR001737",  # SAM-dependent MTase
        "IPR002295",  # Restriction enzyme
        "IPR011991",  # SAM-MTase superfamily
        "IPR025743",  # MTase type 12
        "IPR025942",  # MTase type 13
        "IPR000690",  # Type II restriction enzyme
        "IPR005135",  # DNA specificity domain

        # GO terms
        "GO:0009038",  # restriction enzyme activity
        "GO:0008168",  # methyltransferase activity
        "GO:0009307",  # DNA restriction-modification system
        "GO:0006306"   # DNA methylation
    ]

    # 2) Methyltransferases *OR* DNA-Binding/Regulatory Proteins
    search_terms_mtase_plus_regulators = search_terms_mtase + [
        # Additional regulatory / transcription factor terms
        "regulator of gene expression", "transcription factor", "transcriptional regulator",
        "repressor", "activator", "coactivator", "corepressor",

        # DNA-binding domain types
        "DNA-binding", "helix-turn-helix", "HTH domain", "winged helix domain",
        "zinc finger", "leucine zipper", "basic helix-loop-helix",

        # Common regulatory families & InterPro
        "ArsR", "SmtB", "LuxR", "two-component regulator", "sigma factor",
        "IPR001845",   # HTH ArsR-type
        "IPR036388",   # Winged helix-like domain
        "IPR001867",   # Helix-turn-helix MerR-type
        "IPR036390",   # Another winged helix superfamily
        "IPR036389",   # Helix-turn-helix superfamily

        # GO terms for DNA-binding & regulation
        "GO:0003700",  # DNA-binding transcription factor activity
        "GO:0006351",  # transcription, DNA-templated
        "GO:0006355",  # regulation of transcription, DNA-templated
        "GO:0006357"   # regulation of transcription by RNA polymerase II
    ]

    if args_type == "MTASE":
        return search_terms_mtase
    elif args_type == "ALL":
        return search_terms_mtase_plus_regulators
    else:
        print(f"Error: Invalid --type: {args_type}. Valid options are 'MTASE' or 'ALL'.", file=sys.stderr)
        sys.exit(1)


def process_tsv(input_file, output_file, search_terms):
    """
    Process a TSV file:
      - Load the file
      - Drop 2nd (index=1) and 15th (index=14) columns if present
      - Filter rows by search_terms (case-insensitive)
      - Reformat the first column if it matches 'gnl|extdb|...'
      - Save filtered results to output
    """
    try:
        # Load TSV
        df = pd.read_csv(input_file, sep="\t")

        # Drop 2nd (index=1) and 15th (index=14) columns if present
        cols_to_drop = []
        if len(df.columns) > 1:
            cols_to_drop.append(df.columns[1])
        if len(df.columns) > 14:
            cols_to_drop.append(df.columns[14])
        df.drop(columns=cols_to_drop, inplace=True, errors="ignore")

        # Build a case-insensitive regex pattern with OR logic
        search_pattern = '|'.join(search_terms)

        # Filter rows if any cell in that row matches the pattern (case-insensitive)
        matches = df.apply(
            lambda row: row.astype(str).str.contains(search_pattern, case=False, na=False).any(),
            axis=1
        )
        filtered_df = df[matches]

        # Reformat first column if there's any data
        if not filtered_df.empty:
            filtered_df.iloc[:, 0] = filtered_df.iloc[:, 0].str.extract(
                r'gnl\|extdb\|(.+)', expand=False
            ).fillna(filtered_df.iloc[:, 0])

        # Save to output
        filtered_df.to_csv(output_file, sep="\t", index=False)
        print(f"Filtered results saved to {output_file}")

    except Exception as e:
        print(f"An error occurred while processing {input_file} -> {output_file}:\n{e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Filter a TSV file based on search terms, exclude specific columns, and reformat the first column."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input TSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    parser.add_argument("--type", choices=["ALL", "MTASE"], default="MTASE",
                        help="Choose 'ALL' to filter for methyltransferases + general regulators, "
                             "or 'MTASE' to filter only for methyltransferases (default).")

    args = parser.parse_args()
    search_terms = get_search_terms(args.type)
    process_tsv(args.input, args.output, search_terms)


if __name__ == "__main__":
    main()
