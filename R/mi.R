# Spectral overlap mutual information formulas.

overlap_mi_from_gamma <- function(gamma) {
  gamma_use <- gamma[is.finite(gamma) & gamma > 0]
  if (length(gamma_use) == 0) return(0)
  gamma_use <- pmin(pmax(gamma_use, 0), 1 - 1e-12)
  -0.5 * sum(log(1 - gamma_use))
}

theoretical_gamma <- function(d, c_ratio, side = c("cells", "genes"), r = 10) {
  side <- match.arg(side)
  d_use <- d[is.finite(d) & d > 0]
  d_use <- d_use[seq_len(min(r, length(d_use)))]

  gamma_out <- rep(NA_real_, r)
  if (length(d_use) == 0) return(gamma_out)

  if (side == "cells") {
    g <- (d_use^4 * c_ratio - 1) / (d_use^2 * (c_ratio * d_use^2 + 1))
  } else {
    g <- (d_use^4 * c_ratio - 1) / (c_ratio * d_use^2 * (d_use^2 + 1))
  }

  g[d_use^4 <= 1 / c_ratio] <- 0
  g <- pmin(pmax(g, 0), 1 - 1e-12)
  gamma_out[seq_along(g)] <- g
  gamma_out
}

theoretical_overlap_mi <- function(
  d_cur,
  d_ref,
  c_cur,
  c_ref = c_cur,
  side = c("cells", "genes"),
  r = 10
) {
  side <- match.arg(side)
  gamma_cur <- theoretical_gamma(d_cur, c_cur, side = side, r = r)
  gamma_ref <- theoretical_gamma(d_ref, c_ref, side = side, r = r)

  gamma_cur[!is.finite(gamma_cur)] <- 0
  gamma_ref[!is.finite(gamma_ref)] <- 0
  gamma_double <- pmin(pmax(gamma_cur * gamma_ref, 0), 1 - 1e-12)

  list(
    mi = overlap_mi_from_gamma(gamma_double),
    gamma_cur = gamma_cur,
    gamma_ref = gamma_ref,
    gamma_double = gamma_double
  )
}

theoretical_overlap_mi_from_eigenvalues <- function(
  cur_eigenvalues,
  ref_eigenvalues,
  n_cur,
  p_cur,
  n_ref = n_cur,
  p_ref = p_cur,
  noise_scale_ref = NULL,
  side = c("cells", "genes"),
  r = 10,
  R = 50,
  p_sim = 800,
  seed = 1
) {
  side <- match.arg(side)
  c_ref <- p_ref / n_ref
  if (is.null(noise_scale_ref)) {
    noise_scale_ref <- estimate_noise_scale(
      ref_eigenvalues,
      c_ratio = c_ref,
      R = R,
      p_sim = p_sim,
      seed = seed
    )
  }

  ref_stats <- estimate_spikes(
    ref_eigenvalues,
    n = n_ref,
    p = p_ref,
    noise_scale = noise_scale_ref,
    r_keep = r,
    R = R,
    p_sim = p_sim,
    seed = seed
  )
  cur_stats <- estimate_spikes(
    cur_eigenvalues,
    n = n_cur,
    p = p_cur,
    noise_scale = noise_scale_ref,
    r_keep = r,
    R = R,
    p_sim = p_sim,
    seed = seed
  )

  mi <- theoretical_overlap_mi(
    d_cur = cur_stats$d_spike,
    d_ref = ref_stats$d_spike,
    c_cur = cur_stats$c_ratio,
    c_ref = ref_stats$c_ratio,
    side = side,
    r = r
  )

  list(
    mi = mi$mi,
    gamma_cur = mi$gamma_cur,
    gamma_ref = mi$gamma_ref,
    gamma_double = mi$gamma_double,
    noise_scale_ref = noise_scale_ref,
    cur_stats = cur_stats,
    ref_stats = ref_stats
  )
}

right_singular_vectors <- function(
  x,
  r = 10,
  drop_first = FALSE,
  use_irlba = TRUE,
  maxit = 2000,
  work = NULL
) {
  x <- as_dense_matrix(x)
  r_need <- as.integer(r) + as.integer(drop_first)
  r_eff <- min(r_need, nrow(x) - 1L, ncol(x) - 1L)
  if (!is.finite(r_eff) || r_eff < 1L) return(NULL)

  if (isTRUE(use_irlba) && requireNamespace("irlba", quietly = TRUE)) {
    if (is.null(work)) work <- max(3L * r_eff, 30L)
    fit <- irlba::irlba(x, nv = r_eff, nu = 0, maxit = maxit, work = work)
    zhat <- fit$v[, seq_len(r_eff), drop = FALSE]
  } else {
    fit <- svd(x, nu = 0, nv = r_eff)
    zhat <- fit$v[, seq_len(r_eff), drop = FALSE]
  }
  if (isTRUE(drop_first) && ncol(zhat) >= 2L) {
    zhat <- zhat[, -1L, drop = FALSE]
  }
  zhat <- zhat[, seq_len(min(as.integer(r), ncol(zhat))), drop = FALSE]

  rownames(zhat) <- colnames(x)
  colnames(zhat) <- paste0("z", seq_len(ncol(zhat)))
  zhat
}

label_target_matrix <- function(labels, center = TRUE, drop_empty = TRUE) {
  if (is.null(labels)) stop("labels must not be NULL.", call. = FALSE)
  labels <- factor(labels)
  if (drop_empty) labels <- droplevels(labels)
  if (nlevels(labels) < 2L) {
    stop("labels must contain at least two classes.", call. = FALSE)
  }

  y <- t(stats::model.matrix(~ labels - 1))
  rownames(y) <- sub("^labels", "", rownames(y))
  colnames(y) <- names(labels)
  if (is.null(colnames(y))) colnames(y) <- as.character(seq_along(labels))

  if (isTRUE(center)) {
    y <- y - rowMeans(y)
  }
  y[!is.finite(y)] <- 0
  y
}

subspace_overlap_mi <- function(z1, z2, eps = 1e-12) {
  z1 <- as_dense_matrix(z1)
  z2 <- as_dense_matrix(z2)

  if (!is.null(rownames(z1)) && !is.null(rownames(z2))) {
    common_cells <- intersect(rownames(z1), rownames(z2))
    if (length(common_cells) < 5L) {
      return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L))
    }
    z1 <- z1[common_cells, , drop = FALSE]
    z2 <- z2[common_cells, , drop = FALSE]
  }

  if (nrow(z1) != nrow(z2)) {
    stop("z1 and z2 must have the same number of rows/cells.", call. = FALSE)
  }
  if (ncol(z1) < 1L || ncol(z2) < 1L) {
    return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  z1 <- scale(z1, center = TRUE, scale = FALSE)
  z2 <- scale(z2, center = TRUE, scale = FALSE)
  z1[!is.finite(z1)] <- 0
  z2[!is.finite(z2)] <- 0

  q1 <- qr.Q(qr(z1))
  q2 <- qr.Q(qr(z2))
  r_eff <- min(ncol(q1), ncol(q2))
  if (r_eff < 1L) {
    return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  cmat <- crossprod(q1, q2)
  gamma <- svd(cmat, nu = 0, nv = 0)$d^2
  gamma <- gamma[seq_len(min(r_eff, length(gamma)))]
  gamma <- pmin(pmax(gamma, 0), 1 - eps)
  list(
    mi = -0.5 * sum(log1p(-gamma)),
    gamma = gamma,
    r_eff = length(gamma)
  )
}

target_subspace_mi <- function(zhat, y, r = 10, drop_first = FALSE, eps = 1e-12, use_irlba = TRUE) {
  ztarget <- right_singular_vectors(y, r = r, drop_first = drop_first, use_irlba = use_irlba)
  if (is.null(ztarget)) {
    return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L, ztarget = NULL))
  }
  out <- subspace_overlap_mi(zhat, ztarget, eps = eps)
  out$ztarget <- ztarget
  out
}

label_subspace_mi <- function(zhat, labels, r = 10, eps = 1e-12, use_irlba = TRUE) {
  y <- label_target_matrix(labels, center = TRUE)
  target_subspace_mi(zhat, y, r = r, drop_first = FALSE, eps = eps, use_irlba = use_irlba)
}

matrix_label_subspace_mi <- function(x, labels, r = 10, eps = 1e-12, use_irlba = TRUE) {
  zhat <- right_singular_vectors(x, r = r, drop_first = FALSE, use_irlba = use_irlba)
  if (is.null(zhat)) {
    return(list(mi = NA_real_, gamma = numeric(0), r_eff = 0L, zhat = NULL, ztarget = NULL))
  }
  out <- label_subspace_mi(zhat, labels, r = r, eps = eps, use_irlba = use_irlba)
  out$zhat <- zhat
  out
}
