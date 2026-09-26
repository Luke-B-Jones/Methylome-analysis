# ONT-based Methylomics analysis scripts

This repository contains scripts used in the **Bagby Lab** for microbial methylomics research and methyltransferase (MTase) discovery.

The code was developed during the associated research projects and is provided to document the analyses performed and support reproducibility.

## Important

This repository is **not an end-to-end pipeline**.

The scripts were developed and run **manually, section by section**, depending on the analysis being performed. The main Bash script should therefore be treated as an analysis record and collection of commands rather than a workflow intended to be run from start to finish.

Users wishing to reuse the code should inspect the relevant section, required inputs, paths, and dependencies before running it.

## Analyses represented

The repository includes scripts used for:

- 6mA, 4mC and 5mC methylation analysis;
- differential methylation analysis with `modkit`;
- methylation entropy analysis;
- genomic feature and operon annotation;
- integration with transcriptomic data;
- motif analysis;
- methyltransferase discovery;
- statistical analyses and figure generation.

## 4mC / 5mC separation

The original cytosine-modification data contained combined 4mC and 5mC calls.

Later in the analysis these were separated and relevant downstream analyses were repeated independently:

- **5mC**: modification code `m`
- **4mC**: modification code `21839`

The reconstructed main Bash script reflects this later separation while otherwise retaining the structure of the original analysis.

## Repository use

Some paths, sample names, filenames, and external resources are specific to the original project and computational environment.

The code is provided primarily as a record of the analyses used in the study, rather than as a general-purpose methylomics software package.

