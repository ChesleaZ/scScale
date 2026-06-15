#!/usr/bin/env Rscript

source("scripts/core/empirical_spectrum.R")

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

input <- arg_value("--input")
input_type <- arg_value("--input-type", "matrix")
out_file <- arg_value("--out", "outputs/empirical_spectrum/eigenvalues.csv")
sep <- arg_value("--sep", "auto")
orientation <- arg_value("--orientation", "cells_by_features")
center <- arg_value("--center", "true") %in% c("TRUE", "true", "1", "yes")
scale <- arg_value("--scale", "false") %in% c("TRUE", "true", "1", "yes")
rank_arg <- arg_value("--rank", NA)
rank <- if (is.na(rank_arg)) NULL else as.integer(rank_arg)
reduction <- arg_value("--reduction", "pca")

if (is.null(input) || !nzchar(input)) {
  stop("Usage: compute_empirical_spectrum.R --input=<path> --input-type=<matrix|pca_scores|seurat|eigenvalues> --out=<csv>", call. = FALSE)
}

if (sep == "auto") {
  sep <- if (grepl("\\.tsv$|\\.txt$", input, ignore.case = TRUE)) "\t" else ","
}

if (input_type == "matrix") {
  mat <- read.table(input, sep = sep, header = TRUE, row.names = 1, check.names = FALSE)
  eig <- covariance_eigenvalues(
    mat,
    orientation = orientation,
    center = center,
    scale = scale,
    rank = rank
  )
  eig$source_type <- "matrix"
} else if (input_type == "pca_scores") {
  scores <- read.table(input, sep = sep, header = TRUE, row.names = 1, check.names = FALSE)
  eig <- pca_score_eigenvalues(scores, center = center)
} else if (input_type == "seurat") {
  obj <- readRDS(input)
  eig <- seurat_pca_eigenvalues(obj, reduction = reduction)
} else if (input_type == "eigenvalues") {
  eig <- standardize_eigenvalue_table(read.csv(input, check.names = FALSE))
} else {
  stop("Unsupported --input-type: ", input_type, call. = FALSE)
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(eig, out_file, row.names = FALSE)
write.csv(eigenvalue_summary(eig$eigenvalue), sub("\\.csv$", "_summary.csv", out_file), row.names = FALSE)

cat("Wrote empirical spectrum to", out_file, "\n")
print(utils::head(eig, 10))
