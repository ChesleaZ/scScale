#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

source("packages/scScale/R/mp.R")
source("packages/scScale/R/bbp.R")
source("packages/scScale/R/mi.R")
source("packages/scScale/R/io.R")

eigen_file <- arg_value("--eigenvalues", "outputs/comparisons/gse123025_mp_bbp_server/eigenvalues.csv")
out_dir <- arg_value("--out-dir", "outputs/comparisons/scscale_vs_github")
n <- as.integer(arg_value("--n", "600"))
p <- as.integer(arg_value("--p", "300"))
r <- as.integer(arg_value("--r", "5"))
side <- arg_value("--side", "cells")
mp_R <- as.integer(arg_value("--mp-R", "6"))
mp_p_sim <- as.integer(arg_value("--mp-p-sim", "250"))
seed <- as.integer(arg_value("--seed", "20260611"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Oracle functions copied directly from the GitHub scripts' median-MP/BBP/MI
# section, with only names changed to avoid collisions.
gh_mp_median_sim <- function(c_ratio, R = 50, p_sim = 800, seed = 1) {
  set.seed(seed)
  n_sim <- round(p_sim / c_ratio)
  n_sim <- max(n_sim, 20)

  meds <- replicate(R, {
    G <- matrix(rnorm(p_sim * n_sim), nrow = p_sim, ncol = n_sim)
    ev <- eigen((G %*% t(G)) / n_sim, only.values = TRUE)$values
    median(ev)
  })
  median(meds)
}

gh_calculate_d_from_mu_cov <- function(mu, J, n, p) {
  c <- p / n
  mu_tilde <- mu / J
  term <- mu_tilde - (1 + c)
  disc <- term^2 - 4 * c

  d <- rep(NA_real_, length(mu))
  ok <- is.finite(disc) & disc >= 0
  d[ok] <- sqrt(pmax((term[ok] + sqrt(disc[ok])) / (2 * c), 0))
  d
}

gh_estimate_spikes_from_eigenvalues <- function(mu, J_hat_fix, n, p, r_keep = 10) {
  c_val <- p / n
  mu <- sort(mu[is.finite(mu) & mu > 0], decreasing = TRUE)
  mp_upper <- J_hat_fix * (1 + sqrt(c_val))^2

  d_all <- gh_calculate_d_from_mu_cov(mu, J_hat_fix, n, p)
  d_spike <- d_all
  d_spike[mu <= mp_upper] <- NA_real_

  d_keep <- d_spike[is.finite(d_spike) & d_spike > 0]
  d_keep <- d_keep[1:min(r_keep, length(d_keep))]

  d_keep_pad <- rep(NA_real_, r_keep)
  if (length(d_keep) > 0) d_keep_pad[1:length(d_keep)] <- d_keep

  list(
    c_val = c_val,
    mu = mu,
    mp_upper = mp_upper,
    d_all = d_all,
    d_spike = d_spike,
    n_spikes = sum(mu > mp_upper, na.rm = TRUE),
    d_keep = d_keep,
    d_keep_pad = d_keep_pad
  )
}

gh_theoretical_gamma_single <- function(d, c, side = c("cells", "genes"), r = 10) {
  side <- match.arg(side)
  d_use <- d[is.finite(d) & d > 0]
  d_use <- d_use[1:min(r, length(d_use))]

  gamma_out <- rep(NA_real_, r)
  if (length(d_use) == 0) return(gamma_out)

  if (side == "cells") {
    g <- (d_use^4 * c - 1) / (d_use^2 * (c * d_use^2 + 1))
  } else {
    g <- (d_use^4 * c - 1) / (c * d_use^2 * (d_use^2 + 1))
  }
  g[d_use^4 <= 1 / c] <- 0
  g <- pmin(pmax(g, 0), 1 - 1e-12)

  gamma_out[1:length(g)] <- g
  gamma_out
}

gh_theoretical_overlap <- function(d_cur, d_ref, c_cur, c_ref = c_cur, side = c("cells", "genes"), r = 10) {
  side <- match.arg(side)
  gamma_cur <- gh_theoretical_gamma_single(d_cur, c_cur, side = side, r = r)
  gamma_ref <- gh_theoretical_gamma_single(d_ref, c_ref, side = side, r = r)
  gamma_cur[!is.finite(gamma_cur)] <- 0
  gamma_ref[!is.finite(gamma_ref)] <- 0
  gamma_double <- gamma_cur * gamma_ref
  gamma_double <- pmin(pmax(gamma_double, 0), 1 - 1e-12)
  gamma_use <- gamma_double[is.finite(gamma_double) & gamma_double > 0]
  mi <- if (length(gamma_use) == 0) 0 else -0.5 * sum(log(1 - gamma_use))

  list(mi = mi, gamma_cur = gamma_cur, gamma_ref = gamma_ref, gamma_double = gamma_double)
}

max_abs_diff <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(0)
  max(abs(a[ok] - b[ok]))
}

mu <- read_eigenvalues_csv(eigen_file)
c_ratio <- p / n

gh_mp_med <- gh_mp_median_sim(c_ratio, R = mp_R, p_sim = mp_p_sim, seed = seed)
pkg_mp_med <- mp_median_sim(c_ratio, R = mp_R, p_sim = mp_p_sim, seed = seed)

gh_noise <- median(mu[is.finite(mu) & mu > 0]) / gh_mp_med
pkg_noise <- estimate_noise_scale(mu, c_ratio = c_ratio, R = mp_R, p_sim = mp_p_sim, seed = seed)

gh_spikes <- gh_estimate_spikes_from_eigenvalues(mu, gh_noise, n = n, p = p, r_keep = r)
pkg_spikes <- estimate_spikes(mu, n = n, p = p, noise_scale = pkg_noise, r_keep = r)

gh_mi <- gh_theoretical_overlap(
  gh_spikes$d_spike,
  gh_spikes$d_spike,
  c_cur = c_ratio,
  c_ref = c_ratio,
  side = side,
  r = r
)
pkg_mi <- theoretical_overlap_mi_from_eigenvalues(
  cur_eigenvalues = mu,
  ref_eigenvalues = mu,
  n_cur = n,
  p_cur = p,
  n_ref = n,
  p_ref = p,
  noise_scale_ref = pkg_noise,
  side = side,
  r = r,
  R = mp_R,
  p_sim = mp_p_sim,
  seed = seed
)

summary <- data.frame(
  eigenvalue_file = eigen_file,
  n = n,
  p = p,
  c_ratio = c_ratio,
  side = side,
  r = r,
  mp_R = mp_R,
  mp_p_sim = mp_p_sim,
  seed = seed,
  github_mp_median = gh_mp_med,
  package_mp_median = pkg_mp_med,
  diff_mp_median = abs(gh_mp_med - pkg_mp_med),
  github_noise_scale = gh_noise,
  package_noise_scale = pkg_noise,
  diff_noise_scale = abs(gh_noise - pkg_noise),
  github_mp_upper = gh_spikes$mp_upper,
  package_mp_upper = pkg_spikes$mp_upper,
  diff_mp_upper = abs(gh_spikes$mp_upper - pkg_spikes$mp_upper),
  github_n_spikes = gh_spikes$n_spikes,
  package_n_spikes = pkg_spikes$n_spikes,
  diff_n_spikes = gh_spikes$n_spikes - pkg_spikes$n_spikes,
  max_abs_diff_d_spike = max_abs_diff(gh_spikes$d_spike, pkg_spikes$d_spike),
  max_abs_diff_gamma_cur = max_abs_diff(gh_mi$gamma_cur, pkg_mi$gamma_cur),
  max_abs_diff_gamma_double = max_abs_diff(gh_mi$gamma_double, pkg_mi$gamma_double),
  github_mi = gh_mi$mi,
  package_mi = pkg_mi$mi,
  diff_mi = abs(gh_mi$mi - pkg_mi$mi)
)

utils::write.csv(summary, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(
    rank = seq_along(gh_spikes$mu),
    eigenvalue = gh_spikes$mu,
    github_d = gh_spikes$d_spike,
    package_d = pkg_spikes$d_spike,
    diff_d = gh_spikes$d_spike - pkg_spikes$d_spike
  ),
  file.path(out_dir, "spike_d_comparison.csv"),
  row.names = FALSE
)

print(summary)

if (
  summary$diff_mp_median != 0 ||
    summary$diff_noise_scale != 0 ||
    summary$diff_mp_upper != 0 ||
    summary$diff_n_spikes != 0 ||
    summary$max_abs_diff_d_spike != 0 ||
    summary$max_abs_diff_gamma_cur != 0 ||
    summary$max_abs_diff_gamma_double != 0 ||
    summary$diff_mi != 0
) {
  stop("scScale differs from the copied GitHub formulas; see validation outputs.", call. = FALSE)
}

cat("Exact match to copied GitHub formulas.\n")
