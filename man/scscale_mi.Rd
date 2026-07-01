\name{scscale_mi}
\alias{scscale_mi}
\title{Compare Two Fitted scScale Objects}
\usage{
scscale_mi(fit, target_fit, n_grid = NULL, U_grid = NULL,
  sampling_rates = NULL, combine = TRUE, empirical = TRUE,
  store_empirical_subspaces = FALSE, use_irlba = TRUE, r = NULL,
  P = NULL, eps = 1e-12)
}
\arguments{
\item{fit, target_fit}{\code{scscale_fit} objects.}
\item{n_grid, U_grid, sampling_rates}{Optional cell-number and UMI grids.}
\item{combine}{Whether to return a combined legacy grid when both dimensions are
supplied.}
\item{empirical}{Whether to compute empirical subspace MI from stored matrices.}
\item{store_empirical_subspaces}{Whether to keep empirical subspace matrices.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{r}{Optional number of directions to keep.}
\item{P}{Optional low-rank subspace alignment matrix.}
\item{eps}{Numerical clipping value.}
}
\description{
Compare two fitted objects with \code{\link{scscale_low_rank_mi}} and optional
legacy scaling summaries.
}
\value{
A list containing theory MI and optional empirical or scaling components.
}
\keyword{models}
