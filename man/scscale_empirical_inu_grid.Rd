\name{scscale_empirical_inu_grid}
\alias{scscale_empirical_inu_grid}
\title{Sample an Empirical I(n, U) Grid}
\usage{
scscale_empirical_inu_grid(x, target, n_grid, U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85,
  1.00), r = 10, n_replicates = 10, n_workers = 1, seed = 1,
  input = c("counts", "normalized"), target_input = input,
  target_depth = 1e4, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), center = TRUE, scale = FALSE,
  eps = 1e-12, use_irlba = TRUE)
}
\arguments{
\item{x}{Feature-by-cell source matrix.}
\item{target}{Feature-by-cell target matrix or cell-level target vector.}
\item{n_grid}{Cell numbers to sample.}
\item{U_grid}{Target UMI depths per cell. Provide either \code{U_grid} or
\code{sampling_rates}. If \code{U_grid} is supplied and \code{sampling_rates}
is omitted, the default sampling-rate grid is ignored.}
\item{sampling_rates}{Count downsampling fractions in \code{(0,1]}.}
\item{r}{Fixed rank used for empirical subspace MI.}
\item{n_replicates}{Number of replicate subsamples per grid point.}
\item{n_workers}{Number of parallel workers for empirical grid evaluation.
The default \code{1} evaluates sequentially.}
\item{seed}{Random seed base.}
\item{input, target_input}{Input type for \code{x} and matrix targets.}
\item{target_depth, count_transform, center, scale}{Preprocessing controls passed to empirical MI calculation.}
\item{eps}{Numerical clipping value for MI.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
}
\description{
Sample cells and sequencing depth, then compute fixed-rank empirical subspace
MI for each \code{(n,U)} grid point.  This is the direct empirical counterpart
to spike-model scaling workflows.
}
\value{
A data frame with sampled \code{n}, requested \code{U}, replicate index,
\code{I_empirical}, effective MI rank, and observed UMI summaries.
}
\keyword{models}
