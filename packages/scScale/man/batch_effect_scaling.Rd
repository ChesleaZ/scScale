\name{batch_effect_scaling}
\alias{batch_effect_scaling}
\alias{batch_effect_mi}
\alias{batch_subspace_mi}
\alias{fit_batch_effect_scaling}
\alias{sample_batch_cells}
\alias{summarize_batch_effect_results}
\title{Batch-effect scaling law fitting}
\description{
Estimate biological information from paired cell-side embeddings or paired
feature-by-cell matrices, then fit the two batch-effect scaling laws adapted
from the BMMC CITE-seq validation workflow in
\code{ChesleaZ/Scaling-law-summer-version-}. These functions are part of
\pkg{scScale}.
}
\usage{
batch_subspace_mi(Z1, Z2, eps = 1e-12)

batch_effect_mi(x, y, r_x = 20, r_y = 20,
  transform_x = c("log1p", "pearson", "none"),
  transform_y = c("log1p", "clr", "none"), n_features = NULL,
  min_cells = 10, clip = 10)

sample_batch_cells(meta, batch_col, m_batch, cells_per_batch, seed = NULL,
  batches = NULL)

summarize_batch_effect_results(results)

fit_batch_effect_scaling(summary_df,
  law = c("batch_number", "cells_per_batch"),
  target_col = "mean_I_bio_norm", min_points = 5)
}
\arguments{
\item{Z1, Z2}{Cell-by-component matrices with cell names in row names.}
\item{eps}{Numerical clipping margin for squared canonical correlations.}
\item{x, y}{Feature-by-cell matrices, sparse Matrix objects, or Seurat-style
objects accepted by \code{counts_matrix()}.}
\item{r_x, r_y}{Number of right singular vectors to keep for \code{x} and
\code{y}.}
\item{transform_x, transform_y}{Feature transformations before SVD.}
\item{n_features}{Optional number of highly variable features to keep from
\code{x}.}
\item{min_cells}{Minimum detected cells for feature selection.}
\item{clip}{Pearson residual clipping threshold.}
\item{meta}{Data frame with cell names in row names.}
\item{batch_col}{Column in \code{meta} defining batches.}
\item{m_batch}{Number of batches to sample.}
\item{cells_per_batch}{Cells to sample from each selected batch.}
\item{seed}{Optional sampling seed.}
\item{batches}{Optional fixed batch labels to use.}
\item{results}{Replicate-level data frame containing \code{I_bio} and
\code{I_bio_norm}.}
\item{summary_df}{Summary data frame containing \code{m_batch} or
\code{cells_per_batch} and the target information column.}
\item{law}{Scaling law to fit. \code{"batch_number"} fits
\code{I = I_inf - C * log(1 - A / m)}. \code{"cells_per_batch"} fits
\code{I = I_inf - C * log(1 + A / s)}.}
\item{target_col}{Column of \code{summary_df} to fit.}
\item{min_points}{Minimum finite points required for nonlinear fitting.}
}
\value{
\code{batch_subspace_mi()} and \code{batch_effect_mi()} return lists with
\code{mi}, normalized \code{mi_norm}, squared canonical correlations
\code{gamma}, and effective rank \code{r_eff}. \code{sample_batch_cells()}
returns sampled cells and batches or \code{NULL} when infeasible.
\code{summarize_batch_effect_results()} returns a data frame of replicate
means and standard errors. \code{fit_batch_effect_scaling()} returns a
\code{scsm_batch_scaling} list with fit status, data, predictions,
coefficients, and error metrics.
}
\examples{
set.seed(1)
Z1 <- matrix(rnorm(80), nrow = 20)
Z2 <- Z1 + matrix(rnorm(80, sd = 0.2), nrow = 20)
rownames(Z1) <- rownames(Z2) <- paste0("cell", seq_len(20))
batch_subspace_mi(Z1, Z2)$mi

df <- data.frame(
  m_batch = c(2, 3, 5, 8, 12),
  mean_I_bio_norm = c(0.12, 0.15, 0.18, 0.20, 0.215)
)
fit_batch_effect_scaling(df, law = "batch_number", min_points = 4)$ok
}
