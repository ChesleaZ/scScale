# Public API. The package intentionally uses the GitHub repository's median-MP
# noise fit as the only noise calibration method.

fit_noise <- function(eigenvalues, n = NULL, p = NULL, R = 50, p_sim = 800, seed = 1) {
  n <- eigen_n(eigenvalues, n)
  p <- eigen_p(eigenvalues, p)
  if (is.null(n) || is.null(p)) {
    stop("n and p are required unless eigenvalues came from find_eigenvalues().", call. = FALSE)
  }
  ev <- eigenvalues_vector(eigenvalues)
  c_ratio <- p / n
  scale <- estimate_noise_scale(
    ev,
    c_ratio = c_ratio,
    R = R,
    p_sim = p_sim,
    seed = seed
  )

  structure(
    list(
      method = "median_mp",
      scale = scale,
      n = n,
      p = p,
      c_ratio = c_ratio,
      R = R,
      p_sim = p_sim,
      seed = seed
    ),
    class = "scsm_noise"
  )
}

fit_spikes <- function(
  eigenvalues,
  n = NULL,
  p = NULL,
  noise = NULL,
  r = 10,
  R = 50,
  p_sim = 800,
  seed = 1
) {
  n <- eigen_n(eigenvalues, n)
  p <- eigen_p(eigenvalues, p)
  if (is.null(n) || is.null(p)) {
    stop("n and p are required unless eigenvalues came from find_eigenvalues().", call. = FALSE)
  }
  ev <- eigenvalues_vector(eigenvalues)
  if (is.null(noise)) {
    noise <- fit_noise(ev, n = n, p = p, R = R, p_sim = p_sim, seed = seed)
  }
  noise_scale <- if (inherits(noise, "scsm_noise")) noise$scale else noise

  fit <- estimate_spikes(
    ev,
    n = n,
    p = p,
    noise_scale = noise_scale,
    r_keep = r,
    R = R,
    p_sim = p_sim,
    seed = seed
  )
  fit$noise <- if (inherits(noise, "scsm_noise")) noise else NULL
  fit$r <- r
  class(fit) <- c("scsm_spikes", class(fit))
  fit
}

mi_theory <- function(cur, ref, side = c("cells", "genes"), r = NULL) {
  side <- match.arg(side)
  if (!inherits(cur, "scsm_spikes") || !inherits(ref, "scsm_spikes")) {
    stop("cur and ref must be objects returned by fit_spikes().", call. = FALSE)
  }
  if (is.null(r)) r <- min(cur$r %||% length(cur$d_keep_pad), ref$r %||% length(ref$d_keep_pad))

  out <- theoretical_overlap_mi(
    d_cur = cur$d_spike,
    d_ref = ref$d_spike,
    c_cur = cur$c_ratio,
    c_ref = ref$c_ratio,
    side = side,
    r = r
  )
  out$side <- side
  out$r <- r
  out
}

mi_theory_from_eigenvalues <- function(
  cur_eigenvalues,
  ref_eigenvalues,
  n_cur = NULL,
  p_cur = NULL,
  n_ref = n_cur,
  p_ref = p_cur,
  side = c("cells", "genes"),
  r = 10,
  R = 50,
  p_sim = 800,
  seed = 1
) {
  side <- match.arg(side)
  n_cur <- eigen_n(cur_eigenvalues, n_cur)
  p_cur <- eigen_p(cur_eigenvalues, p_cur)
  n_ref <- eigen_n(ref_eigenvalues, n_ref)
  p_ref <- eigen_p(ref_eigenvalues, p_ref)
  if (is.null(n_cur) || is.null(p_cur) || is.null(n_ref) || is.null(p_ref)) {
    stop("n and p are required unless inputs came from find_eigenvalues().", call. = FALSE)
  }
  cur_ev <- eigenvalues_vector(cur_eigenvalues)
  ref_ev <- eigenvalues_vector(ref_eigenvalues)
  ref_noise <- fit_noise(ref_ev, n = n_ref, p = p_ref, R = R, p_sim = p_sim, seed = seed)
  ref <- fit_spikes(ref_ev, n = n_ref, p = p_ref, noise = ref_noise, r = r)
  cur <- fit_spikes(cur_ev, n = n_cur, p = p_cur, noise = ref_noise, r = r)
  mi <- mi_theory(cur, ref, side = side, r = r)

  list(
    mi = mi$mi,
    gamma_cur = mi$gamma_cur,
    gamma_ref = mi$gamma_ref,
    gamma_double = mi$gamma_double,
    noise = ref_noise,
    cur = cur,
    ref = ref
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
