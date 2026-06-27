\name{scScale-internals}
\alias{scScale-internals}
\alias{discrete_vector_empirical_mi}
\alias{empirical_mi}
\alias{label_target_matrix}
\alias{numeric_target_matrix}
\alias{numeric_vector_empirical_mi}
\alias{right_singular_vectors}
\alias{scscale_cell_scaling}
\alias{scscale_downsample_counts}
\alias{scscale_downsample_counts_fraction}
\alias{scscale_fit_mp_iterative}
\alias{scscale_i_infinity}
\alias{scscale_information_from_theta}
\alias{scscale_mp_goodness}
\alias{scscale_normalize_counts}
\alias{scscale_prepare_matrix}
\alias{scscale_recoverability}
\alias{scscale_spectrum}
\alias{scscale_spike_table}
\alias{scscale_target_mi}
\alias{scscale_theory_mi}
\alias{scscale_theta_infinity}
\alias{scscale_umi_scaling}
\alias{select_hvgs}
\alias{subspace_overlap_mi}
\title{Internal scScale helper functions}
\description{
Internal helper functions used by \code{scscale_fit()}, \code{scscale_mi()},
and the package tutorials. These functions are intentionally not exported.
}
\usage{
discrete_vector_empirical_mi(x, labels, r = 10, eps = 1e-12,
  use_irlba = TRUE)

empirical_mi(x, target, r = 10, eps = 1e-12, use_irlba = TRUE)

label_target_matrix(labels, center = TRUE, drop_empty = TRUE)

numeric_target_matrix(values, center = TRUE, scale = TRUE)

numeric_vector_empirical_mi(x, values, r = 10, eps = 1e-12,
  use_irlba = TRUE)

right_singular_vectors(x, r = 10, drop_first = FALSE, use_irlba = TRUE,
  maxit = 2000, work = NULL)

scscale_cell_scaling(d2_X, p, n_grid, theta_Y, eps = 1e-12)

scscale_downsample_counts(counts, U, seed = NULL)

scscale_downsample_counts_fraction(counts, fraction, seed = NULL)

scscale_fit_mp_iterative(eigenvalues, c_X = NULL, c_ratio = NULL,
  max_iter = 50, grid_n = 20000,
  stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  stop_tol = 1e-4, min_iter = 2)

scscale_i_infinity(x, y = NULL, r = NULL, eps = 1e-12)

scscale_information_from_theta(theta_X, theta_Y, eps = 1e-12)

scscale_mp_goodness(eigenvalues, c_X, tau2, grid_n = 20000,
  eps = 1e-12)

scscale_normalize_counts(counts, target_depth = 1e4, center = TRUE,
  scale = TRUE)

scscale_prepare_matrix(x, input, target_depth = 1e4, center = TRUE,
  scale = TRUE, n_features = NULL, min_cells = 10, n = NULL, p = NULL)

scscale_recoverability(q, c_X = NULL, c_ratio = NULL)

scscale_spectrum(x, input = "normalized", n = NULL, p = NULL)

scscale_spike_table(eigenvalues, bulk, r = NULL)

scscale_target_mi(x, target, r = NULL, eps = 1e-12, use_irlba = TRUE)

scscale_theory_mi(fit, target_fit, r = NULL, eps = 1e-12)

scscale_theta_infinity(q)

scscale_umi_scaling(counts, U_grid = NULL, sampling_rates = NULL, r = 10,
  target_depth = 1e4, n_features = NULL, min_cells = 10, seed = 1,
  n_replicates = 1, theta_Y = NULL, reference_fit = NULL, empirical = TRUE,
  use_irlba = TRUE, mp_max_iter = 50, mp_grid_n = 20000)

select_hvgs(counts, n_features = 2000, min_cells = 10)

subspace_overlap_mi(z1, z2, eps = 1e-12)
}
\arguments{
\item{x, counts}{Feature-by-cell matrix or fitted object.}
\item{labels, values, target}{Empirical target values.}
\item{r}{Number of singular vectors or spikes.}
\item{eps}{Numerical clipping value.}
\item{use_irlba, maxit, work, drop_first}{Singular-vector controls.}
\item{center, scale, target_depth, n_features, min_cells}{Preprocessing controls.}
\item{drop_empty}{Whether to drop empty factor levels.}
\item{d2_X, p, n_grid, theta_Y, theta_X, q}{Spike and MI scaling quantities.}
\item{U, fraction, seed, n_replicates}{Downsampling controls.}
\item{eigenvalues, c_X, c_ratio, tau2, bulk}{MP bulk and spike-fit quantities.}
\item{max_iter, grid_n, stop_metric, stop_tol, min_iter}{MP fitting controls.}
\item{y, fit, target_fit, reference_fit}{Fitted objects or target-side inputs.}
\item{input, n}{Matrix input type and dimensions.}
\item{U_grid, sampling_rates, empirical, mp_max_iter, mp_grid_n}{UMI scaling controls.}
\item{z1, z2}{Cell-by-component subspace matrices.}
}
\value{
Implementation-specific objects used by the exported API.
}
\keyword{internal}
