#!/usr/bin/env Rscript

# Empirical spectral utilities for centered cell-by-feature or feature-by-cell
# matrices. These functions avoid project-specific IO so they can be reused in
# reports, simulations, and Slurm jobs.

as_numeric_matrix <- function(x) {
  if (inherits(x, "Matrix")) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  if (is.matrix(x)) return(x)
  stop("Expected a matrix, sparse Matrix, or data.frame.", call. = FALSE)
}

center_matrix <- function(x, margin = c("columns", "rows")) {
  margin <- match.arg(margin)
  x <- as_numeric_matrix(x)
  if (margin == "columns") {
    if (inherits(x, "sparseMatrix")) {
      means <- Matrix::colMeans(x)
      return(sweep(as.matrix(x), 2, means, "-"))
    }
    return(scale(x, center = TRUE, scale = FALSE))
  }
  if (inherits(x, "sparseMatrix")) {
    means <- Matrix::rowMeans(x)
    return(sweep(as.matrix(x), 1, means, "-"))
  }
  t(scale(t(x), center = TRUE, scale = FALSE))
}

scale_columns <- function(x, center = TRUE, scale = TRUE, scale_floor = 1e-12) {
  x <- as_numeric_matrix(x)
  if (inherits(x, "sparseMatrix")) x <- as.matrix(x)
  if (!center && !scale) return(x)
  out <- if (center) scale(x, center = TRUE, scale = FALSE) else x
  if (scale) {
    s <- apply(out, 2, stats::sd)
    s[!is.finite(s) | s < scale_floor] <- 1
    out <- sweep(out, 2, s, "/")
  }
  out
}

covariance_eigenvalues <- function(
  x,
  orientation = c("cells_by_features", "features_by_cells"),
  center = TRUE,
  scale = FALSE,
  rank = NULL,
  use_irlba = TRUE
) {
  orientation <- match.arg(orientation)
  x <- as_numeric_matrix(x)
  if (orientation == "features_by_cells") x <- t(x)
  if (center || scale) x <- scale_columns(x, center = center, scale = scale)

  n_cells <- nrow(x)
  n_features <- ncol(x)
  denom <- max(1, n_cells - 1)
  max_rank <- min(n_cells - 1L, n_features)
  if (max_rank < 1L) {
    return(data.frame(
      rank = integer(),
      eigenvalue = numeric(),
      n_cells = n_cells,
      n_features = n_features
    ))
  }

  if (is.null(rank) || rank >= max_rank || !use_irlba) {
    sv <- svd(x, nu = 0, nv = 0)$d
  } else {
    if (!requireNamespace("irlba", quietly = TRUE)) {
      stop("Package 'irlba' is required for truncated spectra.", call. = FALSE)
    }
    sv <- irlba::irlba(x, nv = min(rank, max_rank), nu = 0)$d
  }

  eig <- sort((sv^2) / denom, decreasing = TRUE)
  data.frame(
    rank = seq_along(eig),
    eigenvalue = eig,
    log10_eigenvalue = log10(pmax(eig, .Machine$double.xmin)),
    n_cells = n_cells,
    n_features = n_features,
    n_components = length(eig),
    spectrum_is_truncated = length(eig) < max_rank
  )
}

empirical_cdf <- function(x) {
  x <- sort(x[is.finite(x)])
  data.frame(x = x, empirical_cdf = seq_along(x) / length(x))
}

eigenvalue_summary <- function(eigenvalues, thresholds = c(0, 1, 2, 5, 10)) {
  vals <- eigenvalues[is.finite(eigenvalues)]
  positive <- vals[vals > 0]
  out <- data.frame(
    n_eigenvalues = length(vals),
    n_positive = length(positive),
    lambda1 = if (length(vals)) max(vals) else NA_real_,
    median_positive = if (length(positive)) stats::median(positive) else NA_real_,
    mean_positive = if (length(positive)) mean(positive) else NA_real_
  )
  for (thr in thresholds) {
    out[[paste0("n_gt_", gsub("\\.", "_", as.character(thr)))]] <- sum(vals > thr)
  }
  out
}

pca_score_eigenvalues <- function(scores, center = TRUE) {
  scores <- as_numeric_matrix(scores)
  if (center) scores <- scale(scores, center = TRUE, scale = FALSE)
  eig <- apply(scores, 2, stats::var)
  eig <- sort(as.numeric(eig[is.finite(eig) & eig >= 0]), decreasing = TRUE)
  data.frame(
    rank = seq_along(eig),
    eigenvalue = eig,
    log10_eigenvalue = log10(pmax(eig, .Machine$double.xmin)),
    n_cells = nrow(scores),
    n_features = ncol(scores),
    n_components = length(eig),
    spectrum_is_truncated = NA,
    source_type = "pca_scores"
  )
}

pca_stdev_eigenvalues <- function(stdev, n_cells = NA_integer_, n_features = NA_integer_) {
  eig <- sort(as.numeric(stdev)^2, decreasing = TRUE)
  eig <- eig[is.finite(eig) & eig >= 0]
  max_rank <- if (is.finite(n_cells) && is.finite(n_features)) min(n_cells - 1L, n_features) else NA_integer_
  data.frame(
    rank = seq_along(eig),
    eigenvalue = eig,
    log10_eigenvalue = log10(pmax(eig, .Machine$double.xmin)),
    n_cells = n_cells,
    n_features = n_features,
    n_components = length(eig),
    spectrum_is_truncated = if (is.na(max_rank)) NA else length(eig) < max_rank,
    source_type = "pca_stdev"
  )
}

seurat_pca_eigenvalues <- function(obj, reduction = "pca") {
  if (!requireNamespace("SeuratObject", quietly = TRUE) &&
      !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat or SeuratObject is required to read PCA from a Seurat object.", call. = FALSE)
  }
  red <- obj[[reduction]]
  if (is.null(red)) stop("Reduction not found in Seurat object: ", reduction, call. = FALSE)

  stdev <- red@stdev
  embeddings <- tryCatch(SeuratObject::Embeddings(obj, reduction = reduction), error = function(e) NULL)
  loadings <- tryCatch(SeuratObject::Loadings(obj, reduction = reduction), error = function(e) NULL)
  n_cells <- if (is.null(embeddings)) ncol(obj) else nrow(embeddings)
  n_features <- if (is.null(loadings)) NA_integer_ else nrow(loadings)

  out <- pca_stdev_eigenvalues(stdev, n_cells = n_cells, n_features = n_features)
  out$source_type <- "seurat_pca"
  out$reduction <- reduction
  out
}

standardize_eigenvalue_table <- function(df) {
  if (!"eigenvalue" %in% names(df)) {
    candidates <- c("lambda", "eigenvalues", "value")
    hit <- candidates[candidates %in% names(df)][1]
    if (is.na(hit)) stop("Could not find an eigenvalue column.", call. = FALSE)
    df$eigenvalue <- df[[hit]]
  }
  eig <- as.numeric(df$eigenvalue)
  keep <- is.finite(eig) & eig >= 0
  out <- df[keep, , drop = FALSE]
  out$eigenvalue <- eig[keep]
  out <- out[order(out$eigenvalue, decreasing = TRUE), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  out$log10_eigenvalue <- log10(pmax(out$eigenvalue, .Machine$double.xmin))
  if (!"n_cells" %in% names(out)) out$n_cells <- NA_integer_
  if (!"n_features" %in% names(out)) out$n_features <- NA_integer_
  if (!"n_components" %in% names(out)) out$n_components <- nrow(out)
  if (!"spectrum_is_truncated" %in% names(out)) out$spectrum_is_truncated <- NA
  if (!"source_type" %in% names(out)) out$source_type <- "eigenvalues"
  out
}
