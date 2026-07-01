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
  mp_stop_tol = 1e-4, mp_min_iter = 2, r = NULL,
  fit_umi = input[1] == "counts",
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  U_grid = NULL, umi_seed = 1, umi_replicates = 1,
  store_matrix = FALSE, use_irlba = TRUE)
}
\arguments{
\item{x}{Feature-by-cell counts, normalized matrix, representation matrix, or
eigenvalue vector.}
\item{input}{Input type.}
\item{target, target_fit}{Optional target inputs used by legacy MI summaries.}
\item{n, p}{Cell and feature counts, required when \code{input = "eigenvalues"}.}
\item{target_depth, count_transform, center, scale, n_features, min_cells}{Count
preprocessing controls.}
\item{mp_max_iter, mp_grid_n, mp_stop_metric, mp_stop_tol, mp_min_iter}{MAD/MP
bulk and spike-fitting controls.}
\item{r}{Optional number of spikes to keep.}
\item{fit_umi, sampling_rates, U_grid, umi_seed, umi_replicates}{Legacy UMI
scaling controls.}
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
