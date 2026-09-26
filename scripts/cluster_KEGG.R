#!/usr/bin/env Rscript

# Load required libraries and install if necessary
if (!requireNamespace("optparse", quietly = TRUE))
    install.packages("optparse", repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
packages <- c("clusterProfiler", "enrichplot", "ggplot2", "DOSE")
for(pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input TSV file with headers 'KO', 'SA_GeneIDs', and 'log2FoldChange'", metavar="FILE"),
  make_option(c("-o", "--output_prefix"), type="character", default="GSEA_output",
              help="Output file prefix [default %default]", metavar="PREFIX"),
  make_option(c("-c", "--organism"), type="character", default="sau",
              help="KEGG organism code [default %default]", metavar="ORG"),
  make_option(c("-n", "--nperm"), type="integer", default=1000,
              help="Number of permutations for GSEA [default %default]", metavar="NUMBER"),
  make_option(c("--minGSSize"), type="integer", default=10,
              help="Minimum gene set size [default %default]", metavar="NUMBER"),
  make_option(c("--pvalueCutoff"), type="double", default=0.05,
              help="p-value cutoff for GSEA [default %default]", metavar="NUMBER"),
  make_option(c("-t", "--treatment"), type="character", default="Treatment1",
              help="Treatment group name [default %default]", metavar="TREATMENT")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required (--input).", call.=FALSE)
}

# Step 1: Read the input file
cat("Reading input file:", opt$input, "\n")
gsea_data <- read.delim(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Head of imported data:\n")
print(head(gsea_data))
str(gsea_data)

# Check that required columns exist
required_cols <- c("SA_GeneIDs", "log2FoldChange")
if (!all(required_cols %in% colnames(gsea_data))) {
  stop("Input file must have headers 'SA_GeneIDs' and 'log2FoldChange'.", call.=FALSE)
}

# Convert log2FoldChange to numeric and filter out NAs
gsea_data$log2FoldChange <- as.numeric(gsea_data$log2FoldChange)
gsea_data <- gsea_data[!is.na(gsea_data$log2FoldChange), ]

# Create a new column 'SA_Gene' by taking the first SA gene ID from the SA_GeneIDs column
gsea_data$SA_Gene <- sapply(gsea_data$SA_GeneIDs, function(x) {
  sa_ids <- trimws(unlist(strsplit(x, ",")))
  if (length(sa_ids) > 0) {
    return(sa_ids[1])
  } else {
    return(NA)
  }
})

if (all(is.na(gsea_data$SA_Gene))) {
  stop("No valid SA gene IDs could be extracted from the 'SA_GeneIDs' column.", call.=FALSE)
}

# Step 2: Create a ranked gene vector for GSEA using SA gene IDs
gene_vector <- gsea_data$log2FoldChange
names(gene_vector) <- gsea_data$SA_Gene
gene_vector <- sort(gene_vector, decreasing = TRUE)
cat("Ranked gene vector (top 6):\n")
print(head(gene_vector))

# Step 3: Perform GSEA using clusterProfiler
cat("Performing GSEA using organism code:", opt$organism, "\n")
gsea_results <- gseKEGG(geneList     = gene_vector,
                        organism     = opt$organism,
                        nPerm        = opt$nperm,
                        minGSSize    = opt$minGSSize,
                        pvalueCutoff = opt$pvalueCutoff,
                        verbose      = FALSE)

gsea_df <- as.data.frame(gsea_results)
cat("Head of GSEA results:\n")
print(head(gsea_df))

# Order results by adjusted p-value
gsea_df_ordered <- gsea_df[order(gsea_df$p.adjust), ]
cat("Ordered GSEA results (by p.adjust, top 6):\n")
print(head(gsea_df_ordered))

# Simplify pathway descriptions by removing redundant organism information
gsea_df_ordered$SimpleDesc <- gsub(" - Staphylococcus aureus subsp\\. aureus N315 \\(MRSA/VSSA\\)", "", gsea_df_ordered$Description)

# Step 4: Create a dataframe for plotting.
# For each enriched pathway, determine direction based on NES and count the number of core enriched genes.
kegg_df <- gsea_df_ordered
kegg_df$treatment <- opt$treatment
kegg_df$direction <- ifelse(kegg_df$NES > 0, "up", "down")
# Count the number of genes in core enrichment (split by "/")
kegg_df$count <- sapply(kegg_df$core_enrichment, function(x) {
  if (is.na(x)) {
    return(0)
  } else {
    return(length(unlist(strsplit(x, "/"))))
  }
})
# Select only the columns needed for the Python plotting: treatment, KEGG_ID, SimpleDesc, direction, count
plot_df <- kegg_df[, c("treatment", "ID", "SimpleDesc", "direction", "count")]
colnames(plot_df)[colnames(plot_df)=="ID"] <- "KEGG_ID"

# Save the aggregated dataframe to a TSV file for downstream plotting.
output_file <- paste0(opt$output_prefix, "_KEGG_dataframe.tsv")
write.table(plot_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("KEGG pathway dataframe saved to:", output_file, "\n")

# Optionally, you can also keep the dotplot code if desired:
p <- ggplot(kegg_df, aes(x = NES, y = SimpleDesc, size = setSize, color = direction)) +
  geom_point() +
  scale_color_manual(values = c("down" = "red", "up" = "blue")) +
  labs(title = "", x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 10))
dot_pdf <- paste0(opt$output_prefix, "_dotplot.pdf")
ggsave(dot_pdf, plot = p, width = 15, height = 12)
cat("Dotplot of enriched pathways saved to:", dot_pdf, "\n")

cat("GSEA analysis complete.\n")
