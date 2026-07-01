\name{scscale_umi_mi}
\alias{scscale_umi_mi}
\title{Fit UMI-Depth MI Scaling}
\usage{
scscale_umi_mi(pair, x_counts, sampling_rates = c(0.10, 0.20, 0.35,
  0.50, 0.70, 0.85, 1.00), count_transform = pair$count_transform,
  seed = 1, mp_max_iter = 300, mp_grid_n = 2000, use_irlba = TRUE)
}
\arguments{
\item{pair}{A \code{\link{scscale_pair_fit}} object.}
\item{x_counts}{Original X-side count matrix to downsample and refit.}
\item{sampling_rates}{Global UMI sampling fractions in \code{(0, 1]}.}
\item{count_transform}{Count transform used for refits.}
\item{seed}{Random seed for downsampling.}
\item{mp_max_iter, mp_grid_n}{MAD/MP spike-fitting controls for refits.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
}
\description{
Refit the X modality over UMI depth, recompute the X/Y subspace alignment, and
return both the refitted theory MI curve and the refitted \code{I_infinity}
bound.
}
\value{
A \code{scscale_umi_mi} object with \code{curve}, \code{q_by_rate},
\code{fit_by_rate}, and \code{P_by_rate}.
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts[1:100, 1:250]
adt <- gse164378_3p_citeseq_hvg$adt_counts[1:35, 1:250]
pair <- scscale_pair_fit(rna, adt, mp_grid_n = 600, mp_max_iter = 50)
umi <- scscale_umi_mi(pair, rna, sampling_rates = c(0.5, 1),
  mp_grid_n = 500, mp_max_iter = 50)
umi$curve
}
\keyword{models}
