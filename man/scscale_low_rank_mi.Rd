\name{scscale_low_rank_mi}
\alias{scscale_low_rank_mi}
\title{Low-Rank Subspace Mutual Information}
\usage{
scscale_low_rank_mi(theta_X, theta_Y, P = NULL, eps = 1e-12)
}
\arguments{
\item{theta_X, theta_Y}{Recoverability vectors for the X and Y low-rank signal
directions.}
\item{P}{Optional subspace alignment matrix between X and Y. If omitted,
directions are matched by rank.}
\item{eps}{Numerical clipping value.}
}
\description{
Compute the theory MI from X/Y recoverability and optional low-rank subspace
alignment. This is the single theory formula used by the scaling helpers.
}
\value{
A list with \code{mi}, \code{I_theory}, singular values \code{sigma}, squared
canonical correlations \code{gamma}, and effective ranks.
}
\examples{
scscale_low_rank_mi(c(0.5, 0.2), c(0.4, 0.1))$mi
}
\keyword{models}
