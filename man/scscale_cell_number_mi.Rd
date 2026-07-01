\name{scscale_cell_number_mi}
\alias{scscale_cell_number_mi}
\title{Fit Cell-Number MI Scaling}
\usage{
scscale_cell_number_mi(pair, n_grid, q_X = pair$x_fit$spikes$q_X,
  P = pair$P)
}
\arguments{
\item{pair}{A \code{\link{scscale_pair_fit}} object.}
\item{n_grid}{Cell-number grid.}
\item{q_X}{Optional X-side spike-strength vector.}
\item{P}{Optional X/Y low-rank subspace alignment matrix.}
}
\description{
Compute the theory MI curve as the number of cells varies.
}
\value{
A data frame with \code{n}, \code{c_X}, and \code{I_theory}.
}
\keyword{models}
