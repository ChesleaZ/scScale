\name{scscale_empirical_mi}
\alias{scscale_empirical_mi}
\title{Compute Empirical Subspace Mutual Information}
\usage{
scscale_empirical_mi(x, target, r = 10, input = c("counts", "normalized"),
  target_input = input, target_depth = 1e4, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), center = TRUE, scale = FALSE,
  eps = 1e-12, use_irlba = TRUE, store_subspaces = FALSE)
}
\arguments{
\item{x}{Feature-by-cell matrix for the source modality.}
\item{target}{Second feature-by-cell matrix or cell-level target vector.}
\item{r}{Number of singular-vector directions to keep.}
\item{input, target_input}{Input type for \code{x} and matrix targets.}
\item{target_depth, count_transform, center, scale}{Preprocessing controls.}
\item{eps}{Numerical clipping value.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{store_subspaces}{Whether to keep fitted subspace matrices in the output.}
}
\description{
Compute empirical low-rank subspace MI directly from data, either between two
matrices or between one matrix and a cell-level target.
}
\value{
A list containing empirical \code{mi}, canonical-correlation-style
\code{gamma}, and effective rank information.
}
\keyword{models}
