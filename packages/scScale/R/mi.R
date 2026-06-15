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
