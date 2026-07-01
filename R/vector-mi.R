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

scscale_empirical_mi <- function(
  x,
  target,
  r = 10,
  input = c("counts", "normalized"),
  target_input = input,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  eps = 1e-12,
  use_irlba = TRUE,
  store_subspaces = FALSE
) {
  input <- match.arg(input)
  target_input <- match.arg(target_input)
  count_transform <- match.arg(count_transform)
  x <- counts_matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))

  target_is_matrix <- is.matrix(target) || inherits(target, "Matrix")
  if (target_is_matrix) {
    y <- counts_matrix(target)
    if (is.null(colnames(y))) colnames(y) <- paste0("cell_", seq_len(ncol(y)))

    common_cells <- intersect(colnames(x), colnames(y))
    if (length(common_cells) < 2L) {
      if (ncol(x) != ncol(y)) {
        stop("x and target matrix must share cell names or have the same number of columns.", call. = FALSE)
      }
      common_cells <- colnames(x)
      colnames(y) <- common_cells
    }
    x <- x[, common_cells, drop = FALSE]
    y <- y[, common_cells, drop = FALSE]

    X <- if (input == "counts") {
      scscale_normalize_counts(
        x,
        target_depth = target_depth,
        count_transform = count_transform,
        center = center,
        scale = scale
      )
    } else {
      as_dense_matrix(x)
    }
    Y <- if (target_input == "counts") {
      scscale_normalize_counts(
        y,
        target_depth = target_depth,
        count_transform = count_transform,
        center = center,
        scale = scale
      )
    } else {
      as_dense_matrix(y)
    }

    r_use <- min(as.integer(r), ncol(X) - 1L, nrow(X), nrow(Y))
    z_X <- right_singular_vectors(X, r = r_use, use_irlba = use_irlba)
    z_Y <- right_singular_vectors(Y, r = r_use, use_irlba = use_irlba)
    out <- subspace_overlap_mi(z_X, z_Y, eps = eps)
    out$I_empirical <- out$mi
    out$target_type <- "matrix"
    out$n <- ncol(X)
    out$p_X <- nrow(X)
    out$p_Y <- nrow(Y)
    out$r <- r_use
    out$cells <- common_cells
    if (isTRUE(store_subspaces)) {
      out$z_X <- z_X
      out$z_Y <- z_Y
    }
  } else {
    target <- align_vector_to_cells(target, colnames(x), "target")
    X <- if (input == "counts") {
      scscale_normalize_counts(
        x,
        target_depth = target_depth,
        count_transform = count_transform,
        center = center,
        scale = scale
      )
    } else {
      as_dense_matrix(x)
    }
    out <- empirical_mi(X, target, r = r, eps = eps, use_irlba = use_irlba)
    out$I_empirical <- out$mi
    out$n <- ncol(X)
    out$p_X <- nrow(X)
    out$r <- min(as.integer(r), out$r_eff)
    out$cells <- colnames(X)
    if (!isTRUE(store_subspaces)) {
      out$zhat <- NULL
      out$ztarget <- NULL
    }
  }

  class(out) <- c("scscale_empirical_mi", class(out))
  out
}
