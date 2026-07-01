\name{scScale}
\alias{scScale}
\docType{package}
\title{Single-Cell Scaling Law Utilities}
\description{
Gaussian spike fitting and low-rank spectral mutual information utilities for
single-cell scaling-law analyses.
}
\details{
The recommended paired-modality workflow is:
\enumerate{
\item Fit matched modalities once with \code{\link{scscale_pair_fit}}.
\item Fit UMI-depth scaling with \code{\link{scscale_umi_mi}}.
\item Fit cell-number scaling with \code{\link{scscale_cell_number_mi}}.
\item Fit the joint cell-number by UMI grid with
\code{\link{scscale_cell_number_by_umi_mi}}.
\item Compare with direct empirical MI from \code{\link{scscale_empirical_mi}}.
}

All theory curves are routed through \code{\link{scscale_low_rank_mi}}, the
single low-rank subspace-alignment MI formula.
}
\seealso{
\code{\link{scscale_fit}},
\code{\link{scscale_pair_fit}},
\code{\link{scscale_low_rank_mi}},
\code{\link{scscale_umi_mi}},
\code{\link{scscale_cell_number_mi}},
\code{\link{scscale_cell_number_by_umi_mi}},
\code{\link{scscale_empirical_mi}}
}
\keyword{models}
