# Standard eigenvalue intake for single-cell scaling analyses.

find_eigenvalues <- function(
  x,
  n = NULL,
  p = NULL,
  n_features = 2000,
  min_cells = 10,
  transform = c("pearson", "log1p", "none"),
  assay = NULL,
  layer = "counts",
  slot = "counts",
  clip = 10
) {
  transform <- match.arg(transform)

  if (is.numeric(x) && is.null(dim(x))) {
    ev <- clean_eigenvalues(x)
    return(new_eigenvalues(
      eigenvalues = ev,
      n = n,
      p = p %||% length(ev),
      source = "eigenvalues",
      transform = "none",
      features = NULL
    ))
  }

  counts <- counts_matrix(x, assay = assay, layer = layer, slot = slot)
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature_", seq_len(nrow(counts)))
  features <- select_hvgs(counts, n_features = n_features, min_cells = min_cells)
  counts_use <- counts[features, , drop = FALSE]

  X <- switch(
    transform,
    pearson = pearson_residuals(counts_use, clip = clip),
    log1p = scaled_log1p(counts_use),
    none = as_dense_matrix(counts_use)
  )
  ev <- cov_eigenvalues(X)

  new_eigenvalues(
    eigenvalues = ev,
    n = ncol(X),
    p = nrow(X),
    source = class(x)[1],
    transform = transform,
    features = rownames(X)
  )
}

new_eigenvalues <- function(eigenvalues, n, p, source, transform, features) {
  structure(
    list(
      eigenvalues = clean_eigenvalues(eigenvalues),
      n = n,
      p = p,
      source = source,
      transform = transform,
      features = features
    ),
    class = "scsm_eigenvalues"
  )
}

clean_eigenvalues <- function(eigenvalues) {
  sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
}

eigenvalues_vector <- function(x) {
  if (inherits(x, "scsm_eigenvalues")) return(x$eigenvalues)
  clean_eigenvalues(x)
}

eigen_n <- function(x, n = NULL) {
  if (!is.null(n)) return(n)
  if (inherits(x, "scsm_eigenvalues")) return(x$n)
  NULL
}

eigen_p <- function(x, p = NULL) {
  if (!is.null(p)) return(p)
  if (inherits(x, "scsm_eigenvalues")) return(x$p)
  NULL
}
