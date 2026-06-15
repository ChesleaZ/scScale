# Minimal single-cell preprocessing used before spectral analysis.

counts_matrix <- function(x, assay = NULL, layer = "counts", slot = "counts") {
  if (inherits(x, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is required to extract counts from a Seurat object.", call. = FALSE)
    }
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(x)
    return(SeuratObject::GetAssayData(x, assay = assay, layer = layer))
  }

  if (inherits(x, "Assay") || inherits(x, "Assay5")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is required to extract counts from a Seurat assay.", call. = FALSE)
    }
    return(SeuratObject::GetAssayData(x, layer = layer))
  }

  if (is.matrix(x) || inherits(x, "Matrix")) return(x)
  stop("x must be an eigenvalue vector, counts matrix, sparse Matrix, or Seurat object.", call. = FALSE)
}

select_hvgs <- function(counts, n_features = 2000, min_cells = 10) {
  X <- counts
  detected <- row_sums(X > 0)
  log_mean <- row_means_log1p(X)
  log_second <- row_means_log1p_square(X)
  feature_var <- log_second - log_mean^2
  feature_var[!is.finite(feature_var)] <- 0
  feature_var[detected < min_cells] <- 0

  if (is.null(names(feature_var))) names(feature_var) <- rownames(X)
  names(sort(feature_var, decreasing = TRUE))[seq_len(min(n_features, length(feature_var)))]
}

pearson_residuals <- function(counts, clip = 10, pseudo_count = 1e-4) {
  counts_mat <- as_dense_matrix(counts)
  cell_totals <- colSums(counts_mat)
  gene_totals <- rowSums(counts_mat)
  overall_total <- sum(cell_totals)
  if (!is.finite(overall_total) || overall_total <= 0) {
    stop("Count matrix has no positive total count.", call. = FALSE)
  }

  cell_totals[cell_totals == 0] <- 1
  gene_totals[gene_totals == 0] <- 1
  expected <- outer(gene_totals, cell_totals) / overall_total
  residual <- (counts_mat - expected) / sqrt(expected + pseudo_count)
  residual[residual > clip] <- clip
  residual[residual < -clip] <- -clip
  rownames(residual) <- rownames(counts_mat)
  colnames(residual) <- colnames(counts_mat)
  residual
}

scaled_log1p <- function(counts) {
  X <- log1p(as_dense_matrix(counts))
  X <- t(scale(t(X), center = TRUE, scale = TRUE))
  X[!is.finite(X)] <- 0
  X
}
