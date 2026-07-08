\name{scscale_target_fit}
\alias{scscale_target_fit}
\title{Construct a Deterministic Target Subspace Fit}
\usage{
scscale_target_fit(y, input = c("matrix", "labels", "numeric"),
  rank = NULL, center = TRUE, scale = FALSE, drop_empty = TRUE,
  rank_tol = 1e-08, use_irlba = TRUE)
}
\arguments{
\item{y}{Target matrix, label vector, or numeric vector. Matrix inputs are
feature-by-cell.}
\item{input}{Target input type.}
\item{rank}{Optional number of target singular vectors. Defaults to the
numerical rank after preprocessing.}
\item{center, scale}{Whether to center or scale target rows before computing
the target subspace. Label targets are typically centered and not scaled.}
\item{drop_empty}{Whether to drop unused label levels.}
\item{rank_tol}{Numerical rank tolerance.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
}
\description{
Construct a deterministic low-rank target object for use in
\code{\link{scscale_pair_fit}}. Unlike \code{\link{scscale_fit}}, this function
does not fit an MP/BBP spike model. It computes the target right singular
vectors directly and sets target recoverability to one.
}
\value{
A \code{scscale_subspace_fit} object with target matrix \code{matrix}, cell
subspace \code{z}, rank \code{r}, and unit recoverability vectors
\code{theta_X} and \code{theta_infinity}.
}
\examples{
labels <- factor(rep(c("T", "B", "Mono"), length.out = 90))
target <- scscale_target_fit(labels, input = "labels")
target$r
}
\keyword{models}
