#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
library(argparse)

# Define command line argument parser
parser <- ArgumentParser(description = "Extract sparse matrix and metadata from Seurat RDS file's SCT@counts to prepare data for classifier training")
parser$add_argument("--sample",
  type = "character", nargs = "+", default = NULL,
  help = "Specify sample names (e.g., E14.5 E18.5 P3 NC1), if not provided will save all data's SCT@counts"
)
parser$add_argument("--rds_file",
  type = "character", required = TRUE,
  help = "Seurat RDS file path"
)
parser$add_argument("--output_dir",
  type = "character", required = TRUE,
  help = "Output directory path"
)
parser$add_argument("--celltype_col",
  type = "character", default = "celltype_CS",
  help = "Cell type annotation column name, default is celltype_CS"
)

# Parse command line arguments
args <- parser$parse_args()

# Extract parameters
samples <- args$sample
rds_file <- args$rds_file
output_dir <- args$output_dir
celltype_col <- args$celltype_col

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Cannot create output directory: ", output_dir, ", please check path or permissions")
  }
}

# Define log file and start sink
log_file <- file.path(output_dir, "For_Classifier_Preparation.log")
sink(log_file, append = TRUE, split = TRUE)
cat("Output directory:", output_dir, "\n")

# Read RDS file
cat("Reading RDS file:", rds_file, "\n")
seurat_obj <- readRDS(rds_file)

# Check available assays
cat("Available assays:", names(seurat_obj@assays), "\n")

# Extract SCT@counts
cat("Extracting SCT@counts...\n")
if ("SCT" %in% names(seurat_obj@assays)) {
  counts <- seurat_obj@assays$SCT@counts
  cat("SCT@counts dimensions:", dim(counts), "\n")
} else {
  stop("SCT assay does not exist, please check RDS file")
}

# Get metadata
meta <- seurat_obj@meta.data
all_cells <- colnames(counts)
all_genes <- rownames(counts)

# Check if celltype_col exists
if (!(celltype_col %in% colnames(meta))) {
  cat("Warning: specified cell type column ", celltype_col, " does not exist in metadata, classifier labels may not be available\n")
} else {
  cat("Cell type column:", celltype_col, "\n")
  cat("Cell type distribution:\n")
  print(table(meta[[celltype_col]]))
}

# Define sample mapping rules
sample_mapping <- list(
  "E14.5" = list(
    "CS" = c("E14-CS1", "E14-CS2"),
    "WT" = c("E14-WT1", "E14-WT2")
  ),
  "E18.5" = list(
    "CS" = c("E18-CS1", "E18-CS2"),
    "WT" = c("E18-WT1", "E18-WT2")
  ),
  "P3" = list(
    "CS" = c("P3-CS1", "P3-CS2", "P3-CS3"),
    "WT" = c("P3-WT1", "P3-WT2", "P3-WT3")
  ),
  "NC1" = list("E15hiseq", "E17hiseq")
)

# Function: Save sparse matrix and related files, check and create directory before saving
save_sparse_data <- function(counts, genes, cells, meta, prefix, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  mtx_file <- file.path(output_dir, paste0(prefix, "_counts.mtx"))
  cat("Saving sparse matrix to:", mtx_file, "\n")
  writeMM(counts, mtx_file)
  cat("Saving gene names and cell barcodes...\n")
  write.csv(data.frame(Gene = genes), file.path(output_dir, paste0(prefix, "_genes.csv")), row.names = FALSE)
  write.csv(data.frame(Barcode = cells), file.path(output_dir, paste0(prefix, "_barcodes.csv")), row.names = FALSE)
  cat("Saving metadata to:", file.path(output_dir, paste0(prefix, "_meta_data.csv")), "\n")
  write.csv(meta, file.path(output_dir, paste0(prefix, "_meta_data.csv")), row.names = TRUE)
  if (celltype_col %in% colnames(meta)) {
    cat("Saving cell type labels to:", file.path(output_dir, paste0(prefix, "_celltypes.csv")), "\n")
    write.csv(data.frame(Barcode = cells, Celltype = meta[[celltype_col]]),
      file.path(output_dir, paste0(prefix, "_celltypes.csv")),
      row.names = FALSE
    )
  }
}

# Main logic
if (is.null(samples)) {
  # If --sample not provided, save complete SCT@counts
  cat("No sample names provided, saving complete SCT@counts...\n")
  save_sparse_data(counts, all_genes, all_cells, meta, "all_samples", output_dir)
} else {
  # If --sample provided, extract data by sample
  cat("Extracting SCT@counts by sample...\n")
  available_samples <- unique(meta$orig.ident)
  cat("Available samples (orig.ident):", available_samples, "\n")
  for (time_point in samples) {
    if (!time_point %in% names(sample_mapping)) {
      cat("❗️Warning: sample", time_point, "not in mapping rules, skipping\n")
      next