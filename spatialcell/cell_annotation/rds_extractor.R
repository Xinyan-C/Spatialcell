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

# Cell types for analysis
SuScs_celltypes <- c("Immature osteoblast", "Axin2⁺-SuSC", "Pro-inflammatory dural fibroblast", "Chondrocyte","Ectocranial SM", "Pre-osteoblasts", "Outer dural fibroblast",  "Ligament-like mesenchyme","Neural progenitor cell", "Gli1⁺-SuSC", "Capillary endothelial cell","Pericyte", "Mature  osteoblast","Mmp13⁺-Pre-osteoblast",  "Neuron", "Glial cell", "Naive B cell", "M2 Macrophage", "M1 Macrophage", "Proliferating preosteoblast", "Arachnoid fibroblast", "Proliferating B cell", "Proliferating meningeal fibroblast", "Inner dural fibroblast", "Mast cell", "Neutrophil","Osteoclast")

# Function to process cell types
process_cell <- function(counts, meta, celltype_col, cells) {
  read_celltypes <- function() {
    unified_file <- "../utils/color_schemes.py"
    if (file.exists(unified_file)) {
      lines <- readLines(unified_file)
      celltype_line <- grep('c\\("Mast cell"', lines, value = TRUE)
      if (length(celltype_line) > 0) {
        text <- celltype_line[1]
        parts <- unlist(strsplit(text, '"'))
        celltypes <- parts[seq(2, length(parts), by = 2)]
        if (length(celltypes) > 0) {
          return(celltypes)
        }
      }
    }
    return(NULL)
  }
  
  # Check if celltype column exists
  if (!(celltype_col %in% colnames(meta))) {
    cat("Cell type column not found, processing all cells\n")
    return(list(counts = counts, meta = meta, cells = cells))
  }
  
  # Get available cell types
  available_celltypes <- unique(meta[[celltype_col]])
  available_celltypes <- available_celltypes[!is.na(available_celltypes)]
  
  # Check data composition for key cell types
  has_axin2 <- "Axin2⁺-SuSC" %in% available_celltypes
  has_gli1 <- "Gli1⁺-SuSC" %in% available_celltypes
  
  if (has_axin2 && has_gli1) {
    SuSCs <- meta[[celltype_col]] %in% SuScs_celltypes
    cell_indices <- which(SuSCs)
    SuSCs_cells <- cells[cell_indices]
    
    # Extract corresponding counts and metadata
    SuSCs_counts <- counts[, SuSCs_cells]
    SuSCs_meta <- meta[SuSCs_cells, ]
    celltype_counts <- table(SuSCs_meta[[celltype_col]])
    celltype_counts <- celltype_counts[celltype_counts > 0]
    ordered_celltypes <- names(sort(celltype_counts, decreasing = FALSE))
    celltype_names <- read_celltypes()
    
    if (!is.null(celltype_names) && length(celltype_names) == length(ordered_celltypes)) {
      celltype_map <- setNames(celltype_names, ordered_celltypes)    
      SuSCs_meta_map <- SuSCs_meta
      original_celltypes <- as.character(SuSCs_meta[[celltype_col]])
      SuSCs_meta_map[[celltype_col]] <- celltype_map[original_celltypes]
      
      return(list(
        counts = SuSCs_counts, 
        meta = SuSCs_meta_map, 
        cells = SuSCs_cells
      ))
    } else {
      return(list(
        counts = SuSCs_counts, 
        meta = SuSCs_meta, 
        cells = SuSCs_cells
      ))
    }
  } else {
    return(list(counts = counts, meta = meta, cells = cells))
  }
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
  # If --sample not provided, process complete SCT@counts
  cat("No sample names provided, processing complete SCT@counts...\n")
  
  # Apply cell type selection
  processed_data <- process_cell(counts, meta, celltype_col, all_cells)
  
  save_sparse_data(processed_data$counts, all_genes, processed_data$cells, 
                   processed_data$meta, "all_samples", output_dir)
} else {
  # If --sample provided, extract data by sample
  cat("Extracting SCT@counts by sample...\n")
  available_samples <- unique(meta$orig.ident)
  cat("Available samples (orig.ident):", available_samples, "\n")
  
  for (time_point in samples) {
    if (!time_point %in% names(sample_mapping)) {
      cat("Warning: sample", time_point, "not in mapping rules, skipping\n")
      next
    }
    cat("Processing sample:", time_point, "\n")
    
    # Get cells for this time point
    if (is.list(sample_mapping[[time_point]]) && all(c("CS", "WT") %in% names(sample_mapping[[time_point]]))) {
      # Old time points: merge CS and WT
      all_samples <- c(sample_mapping[[time_point]][["CS"]], sample_mapping[[time_point]][["WT"]])
      all_cells_time <- all_cells[meta$orig.ident %in% all_samples]
      sample_desc <- paste("CS/WT samples:", paste(all_samples, collapse = ", "))
    } else if (is.list(sample_mapping[[time_point]])) {
      # NC1: merge orig.ident in list (like E15hiseq, E17hiseq)
      all_samples <- unlist(sample_mapping[[time_point]])
      all_cells_time <- all_cells[meta$orig.ident %in% all_samples]
      sample_desc <- paste("orig.ident samples:", paste(all_samples, collapse = ", "))
    } else {
      # Single orig.ident (reserved for future expansion)
      sample_id <- sample_mapping[[time_point]]
      all_cells_time <- all_cells[meta$orig.ident == sample_id]
      sample_desc <- paste("orig.ident:", sample_id)
    }
    
    if (length(all_cells_time) == 0) {
      cat("Warning: sample", time_point, " (", sample_desc, ") no matching cells found, skipping\n")
    } else {
      cat("Sample", time_point, " (", sample_desc, ") initial matching cells:", length(all_cells_time), "\n")
      
      # Extract subset for this time point
      counts_time <- counts[, all_cells_time]
      meta_time <- meta[all_cells_time, ]
      
      if (celltype_col %in% colnames(meta_time)) {
        cat("Sample", time_point, "initial cell type distribution:\n")
        print(table(meta_time[[celltype_col]]))
      }
      
      # Apply cell type selection for this time point
      cat("Applying cell type selection for sample", time_point, "...\n")
      processed_data <- process_cell(counts_time, meta_time, celltype_col, all_cells_time)
      
      # Save processed data
      save_sparse_data(processed_data$counts, all_genes, processed_data$cells, 
                       processed_data$meta, paste0(time_point, "_all"), output_dir)
    }
  }
}

cat("Complete! Output files saved to:", output_dir, "\n")
# Close sink
sink()