\name{scScale}
\alias{scScale}
\alias{scscale_fit}
\alias{scscale_mi}
\title{Single-Cell Scaling Law Utilities}
\description{
Object-style Gaussian spike fitting and spectral mutual information utilities
for single-cell scaling-law analyses.
}
\usage{
scscale_fit(x, input = c("counts", "normalized", "representation",
  "eigenvalues"), target = NULL, target_fit = NULL, n = NULL, p = NULL,
  target_depth = 1e4, center = TRUE, scale = TRUE, n_features = NULL,
  min_cells = 10, mp_max_iter = 50, mp_grid_n = 20000,
  mp_stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  mp_stop_tol = 1e-4, mp_min_iter = 2, r = NULL,
  fit_umi = input[1] == "counts",
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  U_grid = NULL, umi_seed = 1, umi_replicates = 1,
  store_matrix = FALSE, use_irlba = TRUE)

scscale_mi(fit, target_fit, n_grid = NULL, U_grid = NULL,
  sampling_rates = NULL, combine = TRUE, r = NULL, eps = 1e-12)
}
\arguments{
\item{x}{A feature-by-cell count matrix, sparse \pkg{Matrix}, normalized
feature-by-cell matrix, low-dimensional representation, or covariance
eigenvalue vector.}
\item{input}{Input type. \code{"counts"} applies library-size normalization,
\code{log1p}, and row standardization. \code{"eigenvalues"} requires
\code{n} and \code{p}.}
\item{target}{Optional empirical target.}
\item{target_fit}{Optional target-side \code{scscale_fit} object.}
\item{n, p}{Cell and feature counts, required when fitting directly from
eigenvalues.}
\item{target_depth, center, scale, n_features, min_cells}{Preprocessing
settings for count matrices.}
\item{mp_max_iter, mp_grid_n, mp_stop_metric, mp_stop_tol, mp_min_iter}{MP bulk
fitting controls.}
\item{r}{Number of spikes or singular vectors to keep.}
\item{fit_umi}{Whether \code{scscale_fit()} should also estimate UMI scaling
curves.}
\item{sampling_rates, U_grid}{Sequencing-depth grid as global sampling rates or
expected total UMI counts.}
\item{umi_seed, umi_replicates}{Random seed and replicate count for stochastic
UMI downsampling.}
\item{store_matrix}{Whether to store the normalized matrix inside the returned
fit object.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial singular-vector
computations when available.}
\item{fit}{An \code{scscale_fit} object for the measured modality.}
\item{n_grid}{Optional cell-number grid for MI scaling.}
\item{combine}{If \code{TRUE}, return the full cell-number by UMI grid when
both grids are supplied to \code{scscale_mi()}.}
\item{eps}{Numerical clipping value for MI calculations.}
}
\value{
\code{scscale_fit()} returns an \code{scscale_fit} object with dimensions,
normalization metadata, covariance spectrum, MP bulk fit, spike table, finite
and infinite recoverability vectors, optional UMI scaling, and optional MI
results. \code{scscale_mi()} returns theoretical MI at the fitted point and,
when requested, UMI, cell-number, or full cell-by-UMI scaling tables.
}
\details{
The public workflow is:
\enumerate{
\item Fit each modality with \code{scscale_fit()}.
\item Use one fitted object as the target in \code{scscale_mi()}.
\item Optionally provide cell-number and/or UMI grids to evaluate scaling laws.
}
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts[, 1:200]
adt <- gse164378_3p_citeseq_hvg$adt_counts[, 1:200]

fit_rna <- scscale_fit(rna, r = 3, fit_umi = FALSE, mp_grid_n = 500)
fit_adt <- scscale_fit(adt, r = 3, fit_umi = FALSE, mp_grid_n = 500)
scscale_mi(fit_rna, fit_adt, n_grid = c(100, 200), r = 3)$cell_scaling
}
\keyword{models}
