\name{scScale}
\alias{scScale}
\alias{fit_noise}
\alias{fit_spikes}
\alias{find_eigenvalues}
\alias{mi_theory}
\alias{mi_theory_from_eigenvalues}
\alias{cov_eigenvalues}
\alias{read_counts_csv}
\alias{read_eigenvalues_csv}
\alias{select_hvgs}
\title{Single-cell scaling law utilities}
\description{
Median Marchenko-Pastur noise calibration, BBP spike inversion, and spectral
overlap mutual information utilities for single-cell scaling analyses in
\pkg{scScale}.
The numerical formulas are adapted from
\code{ChesleaZ/Scaling-Laws-in-Representation-Learning}, with an
eigenvalue-first interface. The package intentionally supports only the
median-MP noise fit used by the original repository.
}
\usage{
find_eigenvalues(x, n = NULL, p = NULL, n_features = 2000,
  min_cells = 10, transform = c("pearson", "log1p", "none"),
  assay = NULL, layer = "counts", slot = "counts", clip = 10)

fit_noise(eigenvalues, n, p, R = 50, p_sim = 800, seed = 1)

fit_spikes(eigenvalues, n, p, noise = NULL, r = 10, R = 50,
  p_sim = 800, seed = 1)

mi_theory(cur, ref, side = c("cells", "genes"), r = NULL)

mi_theory_from_eigenvalues(cur_eigenvalues, ref_eigenvalues, n_cur, p_cur,
  n_ref = n_cur, p_ref = p_cur, side = c("cells", "genes"), r = 10,
  R = 50, p_sim = 800, seed = 1)
}
\arguments{
\item{eigenvalues, cur_eigenvalues, ref_eigenvalues}{Covariance eigenvalues,
typically singular values squared divided by cell/sample count. These can also
be objects returned by \code{find_eigenvalues()}.}
\item{x}{A numeric eigenvalue vector, feature-by-cell count matrix, sparse
Matrix, or Seurat object.}
\item{n, n_cur, n_ref}{Number of cells/samples.}
\item{p, p_cur, p_ref}{Number of genes/features.}
\item{n_features}{Number of variable features to use when \code{x} is a
single-cell count object.}
\item{min_cells}{Minimum detected cells for a feature to be eligible.}
\item{transform}{Matrix transform before computing covariance eigenvalues.
\code{"pearson"} follows the original repository's Pearson residual route.}
\item{assay, layer, slot}{Optional Seurat extraction controls.}
\item{clip}{Pearson residual clipping threshold.}
\item{R}{Number of Gaussian simulations used for the median MP estimate.}
\item{p_sim}{Number of simulated features used for the MP median estimate.}
\item{seed}{Random seed for MP median simulation.}
\item{noise}{Optional noise object returned by \code{fit_noise()}, or a numeric
noise scale. If \code{NULL}, median-MP noise is fit from \code{eigenvalues}.}
\item{cur, ref}{Spike-fit objects returned by \code{fit_spikes()}.}
\item{r}{Number of leading spike/overlap dimensions to keep.}
\item{side}{Whether to use the cell-side or gene-side overlap formula.}
}
\value{
Lists or data frames containing the fitted noise scale, MP upper edge, BBP
spike strengths, theoretical overlap gammas, and mutual information.
}
\examples{
mu <- c(12, 6, 3, 1.2, 1, 0.9, 0.8)
noise <- fit_noise(mu, n = 100, p = 50, R = 2, p_sim = 50)
fit <- fit_spikes(mu, n = 100, p = 50, noise = noise)
mi_theory(fit, fit, side = "cells", r = 2)$mi
}
