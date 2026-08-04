\name{scScale-internals}
\alias{scScale-internals}
\alias{discrete_vector_empirical_mi}
\alias{empirical_mi}
\alias{label_target_matrix}
\alias{numeric_target_matrix}
\alias{numeric_vector_empirical_mi}
\alias{right_singular_vectors}
\alias{scscale_downsample_counts}
\alias{scscale_downsample_counts_fraction}
\alias{scscale_fit_mp_iterative}
\alias{scscale_mp_goodness}
\alias{scscale_normalize_counts}
\alias{scscale_prepare_matrix}
\alias{scscale_recoverability}
\alias{scscale_spectrum}
\alias{scscale_spike_table}
\alias{scscale_target_mi}
\alias{scscale_theta_infinity}
\alias{select_hvgs}
\alias{subspace_overlap_mi}
\title{Lower-Level scScale Helper Functions}
\description{
Lower-level helpers used by \code{scscale_fit()}, \code{scscale_mi()}, and the
package tutorials. Most users should start with the functions documented in
\code{\link{scScale}}.
}
\usage{
right_singular_vectors(x, r = 10, drop_first = FALSE, use_irlba = TRUE,
  maxit = 2000, work = NULL)

subspace_overlap_mi(z1, z2, eps = 1e-12)

scscale_downsample_counts(counts, U, seed = NULL)

scscale_downsample_counts_fraction(counts, fraction, seed = NULL)

scscale_fit_mp_iterative(eigenvalues, c_X = NULL, c_ratio = NULL,
  max_iter = 300, grid_n = 20000,
  stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  stop_tol = 1e-4, min_iter = 2,
  median_method = c("sample", "integrate", "grid"),
  mp_sample_p = 1200, mp_sample_reps = 3, mp_sample_seed = 1,
  gap_z_threshold = 6, quiet_z_threshold = 5, quiet_run = 20,
  background_skip = 25, background_window = 500,
  central_quantiles = c(0.10, 0.70))

scscale_mp_goodness(eigenvalues, c_X, tau2, grid_n = 20000,
  eps = 1e-12)

scscale_normalize_counts(counts, target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE, scale = FALSE)

scscale_prepare_matrix(x, input, target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE, scale = FALSE, n_features = NULL, min_cells = 10,
  n = NULL, p = NULL)

scscale_recoverability(q, c_X = NULL, c_ratio = NULL)

scscale_spectrum(x, input = "normalized", n = NULL, p = NULL)

scscale_spike_table(eigenvalues, bulk, r = NULL, force_spikes = FALSE)

scscale_theta_infinity(q)

}
\arguments{
\item{x, counts}{Feature-by-cell matrix or fitted object.}
\item{z1, z2}{Cell-by-component subspace matrices.}
\item{r}{Number of singular vectors or spikes.}
\item{eps}{Numerical clipping value.}
\item{use_irlba, maxit, work, drop_first}{Singular-vector controls.}
\item{center, scale, target_depth, count_transform, n_features, min_cells}{Preprocessing controls.}
\item{q}{Spike and MI scaling quantities.}
\item{U, fraction, seed}{Downsampling controls.}
\item{eigenvalues, c_X, c_ratio, tau2, bulk}{MP bulk and spike-fit quantities.}
\item{force_spikes}{Whether to force the retained top components to be treated
as spikes.}
\item{max_iter, grid_n, stop_metric, stop_tol, min_iter, median_method,
mp_sample_p, mp_sample_reps, mp_sample_seed, gap_z_threshold,
quiet_z_threshold, quiet_run, background_skip, background_window,
central_quantiles}{MP fitting controls.}
\item{input, n}{Matrix input type and dimensions.}
}
\value{
Implementation-specific objects used by the exported API.
}
\keyword{internal}
