\name{scscale_cell_number_by_umi_mi}
\alias{scscale_cell_number_by_umi_mi}
\title{Fit Joint Cell-Number by UMI MI Scaling}
\usage{
scscale_cell_number_by_umi_mi(pair, umi, n_grid,
  sampling_rates = umi$sampling_rates)
}
\arguments{
\item{pair}{A \code{\link{scscale_pair_fit}} object.}
\item{umi}{A \code{\link{scscale_umi_mi}} object.}
\item{n_grid}{Cell-number grid.}
\item{sampling_rates}{Subset or ordering of UMI sampling rates to use.}
}
\description{
Compute the theory MI grid over both cell number and UMI depth, reusing the
refitted UMI spike strengths and subspace alignments.
}
\value{
A data frame with \code{n}, \code{c_X}, \code{I_theory}, and
\code{sampling_rate}.
}
\keyword{models}
