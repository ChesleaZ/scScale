\name{scscale_fit}
\alias{scscale_fit}
\title{Fit a Gaussian Spike Model}
\usage{
scscale_fit(x, input = c("counts", "normalized", "representation",
  "eigenvalues"), target = NULL, target_fit = NULL, n = NULL, p = NULL,
  target_depth = 1e4, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), center = TRUE, scale = FALSE,
  n_features = NULL, min_cells = 10, mp_max_iter = 300, mp_grid_n = 20000,
  mp_stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  mp_stop_tol = 1e-4, mp_min_iter = 2, mp_gap_z_threshold = 6,
  mp_quiet_z_threshold = 5, mp_quiet_run = 20,
  mp_background_skip = 25, mp_background_window = 500,
  mp_central_quantiles = c(0.10, 0.70), r = NULL,
  force_spikes = FALSE,
  fit_umi = input[1] == "counts",
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  U_grid = NULL, umi_linear_intercept = NULL, umi_linear_slope = NULL,
  umi_reference_U = NULL, umi_clamp_rate = TRUE, umi_seed = 1,
  umi_replicates = 1, store_matrix = FALSE, use_irlba = TRUE)
}
\description{
Fit the input matrix with the Gaussian spike model, estimate MP bulk and spike
parameters, and optionally attach target and UMI-scaling summaries.
}
\arguments{
\item{x}{Feature-by-cell counts, normalized matrix, representation matrix, or
eigenvalue vector.}
\item{input}{Input type.}
\item{target, target_fit}{Optional target inputs used by legacy MI summaries.}
\item{n, p}{Cell and feature counts, required when \code{input = "eigenvalues"}.}
\item{target_depth, count_transform, center, scale, n_features, min_cells}{Count
preprocessing controls.}
\item{mp_max_iter, mp_grid_n, mp_stop_metric, mp_stop_tol, mp_min_iter,
mp_gap_z_threshold, mp_quiet_z_threshold, mp_quiet_run, mp_background_skip,
mp_background_window, mp_central_quantiles}{MAD/MP bulk and spike-fitting
controls.}
\item{r}{Optional number of spikes to keep.}
\item{force_spikes}{Whether to force the retained top components to be treated
as spikes.}
\item{fit_umi, sampling_rates, U_grid, umi_seed, umi_replicates}{Legacy UMI
scaling controls.}
\item{umi_linear_intercept, umi_linear_slope, umi_reference_U,
umi_clamp_rate}{Optional controls for linearized UMI-rate calibration.}
\item{store_matrix}{Whether to store the normalized matrix in the returned fit.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
}
\value{
An \code{scscale_fit} object containing the spectrum, MP bulk fit, spike table,
recoverability vectors, and optional stored normalized matrix.
}
\examples{
data(gse123025_myeloid_hvg)
x <- gse123025_myeloid_hvg$counts[1:100, 1:200]
fit <- scscale_fit(x, mp_grid_n = 600, mp_max_iter = 50)
head(fit$spikes)
}
\keyword{models}
