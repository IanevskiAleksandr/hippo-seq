#!/usr/bin/env Rscript

# Load required libraries
required_packages <- c("dplyr", "tidyr", "tibble", "data.table", "openxlsx", 
                       "ggplot2", "patchwork", "tools", "scales", "jsonlite")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Please install it using install.packages('%s')", pkg, pkg))
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

#saveRDS(args, "args.RDS")

if (length(args) != 2) {
  stop("Usage: Rscript CellType.R <input_file> <output_directory>")
}

input_file <- args[1]
output_dir <- args[2]
#input_file = paste0(output_dir, "/", input_file)

# Set output files
markers_file <- file.path("/var/www/html/testhippo", "markersMouseHippo.xlsx")
output_file <- file.path(output_dir, "umap_clusters_and_types.pdf")
json_file <- file.path(output_dir, "results.json")


#' Smart File Reader for Compressed and Uncompressed Files with Command Line Support
#'
#' @description Reads data files from command line arguments, handling both compressed 
#' and uncompressed formats. For zip files, creates a named directory and handles multiple files.
#' 
#' Usage from command line:
#' Rscript script.R input_file.csv
#' Rscript script.R input_file.zip
#' Rscript script.R input_file.gz
#'
#' @return data.frame containing the read data

smart_read <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("File '%s' does not exist", file_path))
  }
  
  tryCatch({
    ext <- tolower(tools::file_ext(file_path))
    base_name <- tools::file_path_sans_ext(basename(file_path))
    
    if (ext %in% c("gz", "zip")) {
      # Create a specific directory for extraction
      extract_dir <- file.path(tempdir(), base_name)
      dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
      on.exit(unlink(extract_dir, recursive = TRUE))
      
      if (ext == "gz") {
        temp_file <- file.path(extract_dir, base_name)
        system_result <- system2("gunzip", 
                                 args = c("-c", file_path), 
                                 stdout = temp_file)
        if (system_result != 0) {
          stop("Failed to decompress gz file")
        }
        return(data.table::fread(temp_file))
      } else { # zip file
        # Unzip to the specific directory
        unzip(file_path, exdir = extract_dir)
        
        # List all files in the extracted directory
        all_files <- list.files(extract_dir, 
                                recursive = TRUE, 
                                full.names = TRUE)
        
        if (length(all_files) == 0) {
          stop("No files found in zip archive")
        }
        
        # Filter for common data file extensions
        data_extensions <- c("csv", "tsv", "txt", "data", "tab")
        data_files <- all_files[tools::file_ext(all_files) %in% data_extensions]
        
        if (length(data_files) == 0) {
          # If no files with common data extensions, use the first file
          selected_file <- all_files[1]
          warning(sprintf("No common data files found. Using first file: %s", 
                          basename(selected_file)))
        } else {
          selected_file <- data_files[1]
          if (length(data_files) > 1) {
            warning(sprintf("Multiple data files found. Using first file: %s", 
                            basename(selected_file)))
          }
        }
        
        message(sprintf("Reading file: %s", basename(selected_file)))
        return(data.table::fread(selected_file))
      }
    }
    
    # For uncompressed files, just use fread directly
    message(sprintf("Reading uncompressed file: %s", basename(file_path)))
    return(data.table::fread(file_path))
    
  }, error = function(e) {
    stop(sprintf("Error reading file '%s': %s", gsub("testhippo","",gsub("html","",gsub("var","",gsub("www","",e$message))))))
  })
}


#' Process Gene Set Data
#'
#' @param markers_file Path to Excel file containing marker genes
#' @param expression_matrix Expression matrix with genes as rows
#' @return List containing processed gene sets and sensitivity scores
process_gene_sets <- function(markers_file, expression_matrix) {
  # Read and process marker genes
  finallistMerged <- openxlsx::read.xlsx(markers_file)
  if (nrow(finallistMerged) == 0) {
    stop("Marker gene file is empty")
  }
  
  tidy_finallist <- finallistMerged %>%
    separate_rows(genes, sep = ",") %>%
    mutate(genes = trimws(genes)) %>%
    as.data.frame()
  
  # Create gene sets
  gs_positive <- split(tidy_finallist$genes, tidy_finallist$celltype)
  
  # Calculate marker sensitivity
  marker_stat <- sort(table(unlist(gs_positive)), decreasing = TRUE)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(
      as.numeric(marker_stat),
      to = c(0, 1),
      from = c(length(gs_positive), 1)
    ),
    gene_ = names(marker_stat),
    stringsAsFactors = FALSE
  )
  
  # Filter expression matrix for marker genes
  all_markers <- unique(tidy_finallist$genes)
  valid_genes <- rownames(expression_matrix) %in% all_markers
  filtered_matrix <- expression_matrix[valid_genes, , drop = FALSE]
  
  # Update gene sets to only include genes present in the expression matrix
  gs_positive <- lapply(gs_positive, function(genes) {
    genes[genes %in% rownames(filtered_matrix)]
  })
  
  # Remove empty gene sets
  gs_positive <- gs_positive[sapply(gs_positive, length) > 0]
  
  if (length(gs_positive) == 0) {
    stop("No valid gene sets remaining after filtering")
  }
  
  return(list(
    gene_sets = gs_positive,
    sensitivity = marker_sensitivity,
    filtered_matrix = filtered_matrix
  ))
}

#' Calculate Enrichment Scores
#'
#' @param Z Scaled expression matrix
#' @param gs List of gene sets
#' @return Matrix of enrichment scores
optimize_scores <- function(Z, gs) {
  if (length(gs) == 0) stop("Empty gene sets provided")
  if (nrow(Z) == 0) stop("Empty expression matrix provided")
  
  es <- matrix(
    nrow = length(gs),
    ncol = ncol(Z),
    dimnames = list(names(gs), colnames(Z))
  )
  
  for (i in seq_along(gs)) {
    gene_indices <- match(gs[[i]], rownames(Z))
    valid_indices <- !is.na(gene_indices)
    
    if (sum(valid_indices) == 0) {
      warning(sprintf("No valid genes found for set %s", names(gs)[i]))
      next
    }
    
    gs_data <- Z[gene_indices[valid_indices], , drop = FALSE]
    es[i, ] <- colSums(gs_data) / sqrt(nrow(gs_data))
  }
  
  return(es)
}

#' Generate UMAP Visualizations
#'
#' @param meta_data Metadata containing UMAP coordinates and clusters
#' @param sctype_scores Cell type assignment scores
#' @param output_file Output file path for the plot
#' @return ggplot object of the combined visualization
generate_umap_plots <- function(meta_data, sctype_scores, output_file = NULL) {
  # Merge metadata with cell type assignments
  meta_data_with_types <- merge(
    meta_data,
    sctype_scores[, c("cluster", "type")],
    by = "cluster",
    all.x = TRUE
  )
  
  # Create cluster plot
  p1 <- ggplot(meta_data_with_types, aes(x = umap_1, y = umap_2, color = factor(cluster))) +
    geom_point(size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Clusters",
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Cluster"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Create cell type plot
  p2 <- ggplot(meta_data_with_types, aes(x = umap_1, y = umap_2, color = type)) +
    geom_point(size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Cell Types",
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Cell Type"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Combine plots
  combined_plot <- p1 + p2
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot, width = 15, height = 6)
    ggsave(gsub("pdf","png",output_file), combined_plot, width = 15, height = 6)
  }
  
  return(combined_plot)
}

#' Main Analysis Pipeline
#'
#' @param expression_file Path to expression data file
#' @param markers_file Path to marker genes file
#' @param output_file Path for output plot
#' @return List containing results and plots
main_analysis <- function(expression_file, markers_file, output_file = "umap_clusters_and_types.pdf") {
  # Read and process expression data
  message("Reading expression data...")
  matr_ <- smart_read(expression_file) %>%
    as.data.frame() %>%
    column_to_rownames(colnames(.)[1])
  
  # Extract metadata
  meta_rows <- which(rownames(matr_) %in% c("cluster", "umap_1", "umap_2"))
  meta_data <- as.data.frame(t(matr_[meta_rows, , drop = FALSE]))
  expr_matrix <- matr_[-meta_rows, ]
  
  # Process gene sets
  message("Processing gene sets...")
  gene_set_data <- process_gene_sets(markers_file, expr_matrix)
  
  # Scale expression data
  message("Scaling expression data...")
  scRNAseqData <- t(scale(t(gene_set_data$filtered_matrix)))
  
  # Get valid genes present in both sensitivity scores and expression data
  valid_genes <- intersect(
    gene_set_data$sensitivity$gene_,
    rownames(scRNAseqData)
  )
  
  if (length(valid_genes) == 0) {
    stop("No overlap between sensitivity scores and expression data genes")
  }
  
  # Filter sensitivity scores for valid genes
  valid_sensitivity <- gene_set_data$sensitivity[
    gene_set_data$sensitivity$gene_ %in% valid_genes,
  ]
  
  # Apply sensitivity scores to valid genes
  Z <- scRNAseqData
  gene_indices <- match(valid_sensitivity$gene_, rownames(Z))
  valid_indices <- !is.na(gene_indices)
  Z[gene_indices[valid_indices], ] <- Z[gene_indices[valid_indices], ] * 
    valid_sensitivity$score_marker_sensitivity[valid_indices]
  
  # Update gene sets to only include genes present in Z
  gs <- lapply(gene_set_data$gene_sets, function(genes) {
    intersect(genes, rownames(Z))
  })
  gs <- gs[sapply(gs, length) > 0]
  
  if (length(gs) == 0) {
    stop("No valid gene sets remaining after filtering")
  }
  
  # Calculate enrichment scores
  message("Calculating enrichment scores...")
  es <- optimize_scores(Z, gs)
  
  # Process cluster results
  cL_results <- do.call("rbind", lapply(unique(meta_data$cluster), function(cl) {
    es.max.cl <- sort(rowSums(es[, rownames(meta_data[meta_data$cluster == cl, ])]), 
                      decreasing = TRUE)
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = es.max.cl,
      ncells = sum(meta_data$cluster == cl)
    ), 10)
  }))
  
  # Get top scores per cluster
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  
  # Mark low-confidence assignments
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
  
  # Generate visualizations
  message("Generating UMAP visualizations...")
  plots <- generate_umap_plots(meta_data, sctype_scores, output_file)
  
  message("Analysis complete!")
  
  return(list(
    sctype_scores = sctype_scores,
    plots = plots,
    enrichment_scores = es,
    gene_sets = gs
  ))
}


# Main execution
tryCatch({
    message(sprintf("Starting analysis with expression file: %s", basename(input_file)))
    
    results <- main_analysis(
      expression_file = input_file,
      markers_file = markers_file,
      output_file = output_file
    )
    
    # Save results to JSON
    sctype_scores_df <- as.data.frame(results$sctype_scores)
    write_json(sctype_scores_df, json_file)
    
    message("Analysis completed successfully!")
    message(sprintf("Results saved in directory: %s", output_dir))
    
}, error = function(e) {
    message(sprintf("Error in execution: %s", e$message))
    quit(status = 1)
})