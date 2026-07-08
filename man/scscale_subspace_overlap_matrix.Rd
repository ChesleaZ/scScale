\name{scscale_subspace_overlap_matrix}
\alias{scscale_subspace_overlap_matrix}
\title{Compute a Subspace Overlap Matrix}
\usage{
scscale_subspace_overlap_matrix(z_X, z_Y)
}
\arguments{
\item{z_X, z_Y}{Cell-by-rank orthonormal subspace bases. If row names are
present, bases are aligned by shared cell names.}
}
\description{
Compute the cross-product matrix between two low-rank cell subspaces after
optional cell-name alignment. This matrix is used by
\code{\link{scscale_low_rank_mi}}.
}
\value{
A numeric overlap matrix with rows corresponding to \code{z_X} directions and
columns corresponding to \code{z_Y} directions.
}
\keyword{models}
