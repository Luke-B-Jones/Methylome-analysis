#!/usr/bin/env Rscript

# Load necessary library
library(DESeq2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# MODIFIED: Now we require at least three arguments: one or more count files, a condition string, and an output file name.
if (length(args) < 3) {
  stop("Usage: script.R <count_file1> [<count_file2> ...] <condition_string> <output_file>")
}

# MODIFIED: The output file name is now provided as the last argument.
output_file <- args[length(args)]
# MODIFIED: The condition string is now the second last argument.
condition_string <- args[length(args) - 1]
# MODIFIED: Count files are all remaining arguments.
count_files <- args[1:(length(args) - 2)]

# Function to load FeatureCounts data and retain only the count column
load_counts <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", comment.char = "#")
  counts <- df[, c("Geneid", names(df)[ncol(df)])]  # Only keep Geneid and final column (counts)
  colnames(counts) <- c("Geneid", basename(file))     # Rename column to file basename for clarity
  return(counts)
}

# Load and merge count data
count_data_list <- lapply(count_files, load_counts)
count_data <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_data_list)
row.names(count_data) <- count_data$Geneid
count_data <- count_data[, -1]  # Remove Geneid column after setting row names

# Parse conditions from input string (e.g., "control,control,control,treatment,treatment,treatment")
conditions <- unlist(strsplit(condition_string, split = ","))
if (length(conditions) != ncol(count_data)) {
  stop("The number of conditions must match the number of columns in the merged count data.")
}

# Prepare colData for DESeq2
col_data <- data.frame(condition = factor(conditions), row.names = colnames(count_data))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Output results to the provided file name
write.csv(as.data.frame(res), file = output_file)

cat("DESeq2 analysis complete. Results saved to", output_file, "\n")
