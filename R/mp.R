# Median Marchenko-Pastur calibration adapted from
# ChesleaZ/Scaling-Laws-in-Representation-Learning.

mp_cache <- new.env(parent = emptyenv())

mp_median_sim <- function(c_ratio, R = 50, p_sim = 800, seed = 1) {
  if (!is.finite(c_ratio) || c_ratio <= 0) {
    stop("c_ratio must be positive.", call. = FALSE)
  }

  set.seed(seed)
  n_sim <- round(p_sim / c_ratio)
  n_sim <- max(n_sim, 20)

  meds <- replicate(R, {
    G <- matrix(stats::rnorm(p_sim * n_sim), nrow = p_sim, ncol = n_sim)
    ev <- eigen((G %*% t(G)) / n_sim, only.values = TRUE)$values
    stats::median(ev)
  })
  stats::median(meds)
}

get_mp_median_cached <- function(c_ratio, R = 50, p_sim = 800, seed = 1) {
  key <- sprintf("c%.6f_R%d_p%d_seed%d", c_ratio, R, p_sim, seed)
  if (exists(key, envir = mp_cache, inherits = FALSE)) {
    return(get(key, envir = mp_cache, inherits = FALSE))
  }
  val <- mp_median_sim(c_ratio = c_ratio, R = R, p_sim = p_sim, seed = seed)
  assign(key, val, envir = mp_cache)
  val
}

estimate_noise_scale <- function(eigenvalues, c_ratio, R = 50, p_sim = 800, seed = 1) {
  mu <- eigenvalues[is.finite(eigenvalues) & eigenvalues > 0]
  if (length(mu) < 2) {
    stop("Need at least two positive covariance eigenvalues.", call. = FALSE)
  }

  med_obs <- stats::median(mu)
  med_theory <- get_mp_median_cached(c_ratio = c_ratio, R = R, p_sim = p_sim, seed = seed)
  med_obs / med_theory
}

estimate_noise_scale_from_matrix <- function(X, R = 50, p_sim = 800, seed = 1) {
  n_real <- ncol(X)
  p_real <- nrow(X)
  estimate_noise_scale(
    cov_eigenvalues(X),
    c_ratio = p_real / n_real,
    R = R,
    p_sim = p_sim,
    seed = seed
  )
}
