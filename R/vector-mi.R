# MI-like task-difficulty scores for vector targets.

align_vector_to_cells <- function(target, cells, target_name = "target") {
  if (is.null(target)) stop(target_name, " must not be NULL.", call. = FALSE)
  if (is.null(cells)) {
    if (length(target) < 1L) stop(target_name, " must be non-empty.", call. = FALSE)
    return(target)
  }

  if (!is.null(names(target))) {
    missing <- setdiff(cells, names(target))
    if (length(missing)) {
      stop(target_name, " is missing values for ", length(missing), " cells.", call. = FALSE)
    }
    return(target[cells])
  }

  if (length(target) != length(cells)) {
    stop(target_name, " length must match the number of cells when names are absent.", call. = FALSE)
  }
  names(target) <- cells
  target
}

numeric_target_matrix <- function(values, center = TRUE, scale = TRUE) {
  values <- as.numeric(values)
  if (anyNA(values)) stop("numeric target contains NA values.", call. = FALSE)
  y <- matrix(values, nrow = 1L)
  rownames(y) <- "numeric_target"
  colnames(y) <- names(values)

  if (isTRUE(center)) y <- y - rowMeans(y)
  if (isTRUE(scale)) {
    sd_y <- stats::sd(as.numeric(y))
    if (!is.finite(sd_y) || sd_y <= 0) {
      stop("numeric target must have positive variance.", call. = FALSE)
    }
    y <- y / sd_y
  }
  y[!is.finite(y)] <- 0
  y
}

discrete_vector_empirical_mi <- function(x, labels, r = 10, eps = 1e-12, use_irlba = TRUE) {
  labels <- align_vector_to_cells(labels, colnames(x), "labels")
  out <- matrix_label_subspace_mi(x, labels, r = r, eps = eps, use_irlba = use_irlba)
  out$target_type <- "discrete"
  out$class_counts <- table(factor(labels))
  out
}

numeric_vector_empirical_mi <- function(x, values, r = 10, eps = 1e-12, use_irlba = TRUE) {
  values <- align_vector_to_cells(values, colnames(x), "values")
  zhat <- right_singular_vectors(x, r = r, drop_first = FALSE, use_irlba = use_irlba)
  if (is.null(zhat)) {
    return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L, zhat = NULL, ztarget = NULL, target_type = "numeric"))
  }
  y <- numeric_target_matrix(values)
  ztarget <- t(y)
  colnames(ztarget) <- "numeric_target"
  out <- subspace_overlap_mi(zhat, ztarget, eps = eps)
  out$zhat <- zhat
  out$ztarget <- ztarget
  out$target_type <- "numeric"
  out$target_mean <- mean(as.numeric(values))
  out$target_sd <- stats::sd(as.numeric(values))
  out
}

empirical_mi <- function(x, target, r = 10, eps = 1e-12, use_irlba = TRUE) {
  if (is.numeric(target) || is.integer(target)) {
    numeric_vector_empirical_mi(x, target, r = r, eps = eps, use_irlba = use_irlba)
  } else {
    discrete_vector_empirical_mi(x, target, r = r, eps = eps, use_irlba = use_irlba)
  }
}
