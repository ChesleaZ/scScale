\name{scaling_fit}
\alias{scaling_fit}
\alias{fit_cell_scaling}
\alias{fit_umi_scaling}
\alias{fit_batch_scaling}
\alias{predict.scscale_fit}
\alias{coef.scscale_fit}
\alias{summary.scscale_fit}
\alias{plot.scscale_fit}
\title{Fit and extrapolate scScale scaling laws}
\description{
Fit R-style scaling-law objects for cell-number, UMI-depth, and batch scaling.
The returned objects share a common \code{scscale_fit} interface supporting
\code{predict()}, \code{coef()}, \code{summary()}, and \code{plot()}.
}
\usage{
fit_cell_scaling(counts, cell_grid, reps = 1, n_features = 300,
  transform = c("log1p", "pearson", "none"), min_cells = 10, r = 8,
  R = 3, p_sim = 100, seed = 1)

fit_umi_scaling(counts, umi_grid, n_cells = min(3000L, ncol(counts_matrix(counts))),
  reps = 1, n_features = 300, transform = c("log1p", "pearson", "none"),
  min_cells = 10, r = 8, R = 3, p_sim = 100, seed = 1)

fit_batch_scaling(summary_df, law = c("batch_number", "cells_per_batch"),
  target_col = "mean_I_bio_norm", min_points = 5)

\method{predict}{scscale_fit}(object, newdata = NULL, ...)

\method{coef}{scscale_fit}(object, ...)

\method{summary}{scscale_fit}(object, ...)

\method{plot}{scscale_fit}(x, xlab = x$x_col, ylab = x$y_col, ...)
}
\arguments{
\item{counts}{Feature-by-cell count matrix, sparse Matrix, or Seurat-style
object accepted by \code{find_eigenvalues()}.}
\item{cell_grid}{Vector of cell counts to subsample.}
\item{umi_grid}{Vector of UMI thinning fractions in \code{(0, 1]}.}
\item{n_cells}{Number of cells to use for UMI-depth fitting.}
\item{reps}{Number of replicate subsamples per grid point.}
\item{n_features}{Number of highly variable features to use.}
\item{transform}{Preprocessing transform passed to \code{find_eigenvalues()}.}
\item{min_cells}{Minimum detected cells for feature selection.}
\item{r}{Retained rank used for spike and MI calculations.}
\item{R, p_sim, seed}{Median-MP simulation and random seed controls.}
\item{summary_df}{Batch-effect summary table, usually from
\code{summarize_batch_effect_results()}.}
\item{law}{Batch law to fit: \code{"batch_number"} or
\code{"cells_per_batch"}.}
\item{target_col}{Column in \code{summary_df} containing the MI response.}
\item{min_points}{Minimum finite points required for nonlinear fitting.}
\item{object, x}{An object returned by one of the fit functions.}
\item{xlab, ylab}{Axis labels for \code{plot()}.}
\item{newdata}{Data frame containing the relevant predictor column, such as
\code{n_cells}, \code{mean_umi_per_cell}, \code{m_batch}, or
\code{cells_per_batch}.}
\item{...}{Additional arguments passed to plotting methods.}
}
\value{
A fit object with class \code{scscale_fit} and a subtype class such as
\code{scscale_cell_fit}, \code{scscale_umi_fit}, or \code{scscale_batch_fit}.
The object stores fitted coefficients, observed grid values, fitted grid
predictions, reference spectral fits, settings, and diagnostics.
}
\examples{
set.seed(1)
counts <- matrix(rpois(60 * 50, lambda = 5), nrow = 60)
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

fit <- fit_cell_scaling(
  counts,
  cell_grid = c(15, 25, 35, 50),
  n_features = 30,
  r = 3,
  R = 2,
  p_sim = 30
)
coef(fit)
predict(fit, data.frame(n_cells = c(60, 80)))
}
