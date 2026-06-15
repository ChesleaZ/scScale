# BBP spike inversion for covariance eigenvalues.

cov_eigenvalues <- function(X) {
  n_real <- ncol(X)
  p_real <- nrow(X)
  sv <- svd(t(X), nu = 0, nv = 0)$d
  (sv[seq_len(min(p_real, n_real))]^2) / n_real
}

bbp_d_from_cov_eigenvalues <- function(eigenvalues, noise_scale, n, p) {
  c_ratio <- p / n
  mu_tilde <- eigenvalues / noise_scale
  term <- mu_tilde - (1 + c_ratio)
  disc <- term^2 - 4 * c_ratio

  d <- rep(NA_real_, length(eigenvalues))
  ok <- is.finite(disc) & disc >= 0
  d[ok] <- sqrt(pmax((term[ok] + sqrt(disc[ok])) / (2 * c_ratio), 0))
  d
}

estimate_spikes <- function(
  eigenvalues,
  n,
  p,
  noise_scale = NULL,
  r_keep = 10,
  R = 50,
  p_sim = 800,
  seed = 1
) {
  mu <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
  if (length(mu) < 2) {
    stop("Need at least two positive covariance eigenvalues.", call. = FALSE)
  }

  c_ratio <- p / n
  if (is.null(noise_scale)) {
    noise_scale <- estimate_noise_scale(mu, c_ratio = c_ratio, R = R, p_sim = p_sim, seed = seed)
  }

  mp_upper <- noise_scale * (1 + sqrt(c_ratio))^2
  d_all <- bbp_d_from_cov_eigenvalues(mu, noise_scale = noise_scale, n = n, p = p)
  d_spike <- d_all
  d_spike[mu <= mp_upper] <- NA_real_

  d_keep <- d_spike[is.finite(d_spike) & d_spike > 0]
  d_keep <- d_keep[seq_len(min(r_keep, length(d_keep)))]

  d_keep_pad <- rep(NA_real_, r_keep)
  if (length(d_keep) > 0) d_keep_pad[seq_along(d_keep)] <- d_keep

  spike_table <- data.frame(
    spike_rank = seq_along(mu),
    eigenvalue = mu,
    noise_scale = noise_scale,
    c_ratio = c_ratio,
    mp_upper = mp_upper,
    is_spike = mu > mp_upper,
    d = d_spike,
    d2 = d_spike^2
  )

  list(
    n = n,
    p = p,
    c_ratio = c_ratio,
    noise_scale = noise_scale,
    mp_upper = mp_upper,
    eigenvalues = mu,
    d_all = d_all,
    d_spike = d_spike,
    n_spikes = sum(mu > mp_upper, na.rm = TRUE),
    d_keep = d_keep,
    d_keep_pad = d_keep_pad,
    q_keep_pad = d_keep_pad^2,
    d1_sq = if (length(d_keep) == 0) NA_real_ else d_keep[1]^2,
    sum_d2 = sum(d_keep^2, na.rm = TRUE),
    sum_d4 = sum(d_keep^4, na.rm = TRUE),
    spike_table = spike_table
  )
}

estimate_spikes_from_matrix <- function(
  X,
  noise_scale = NULL,
  r_keep = 10,
  R = 50,
  p_sim = 800,
  seed = 1
) {
  estimate_spikes(
    cov_eigenvalues(X),
    n = ncol(X),
    p = nrow(X),
    noise_scale = noise_scale,
    r_keep = r_keep,
    R = R,
    p_sim = p_sim,
    seed = seed
  )
}
