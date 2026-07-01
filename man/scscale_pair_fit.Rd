\name{scscale_pair_fit}
\alias{scscale_pair_fit}
\title{Fit Matched Modalities for MI Scaling}
\usage{
scscale_pair_fit(x, y, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), mp_max_iter = 300, mp_grid_n = 3000,
  use_irlba = TRUE, ...)
}
\arguments{
\item{x, y}{Matched feature-by-cell count matrices for the two modalities.}
\item{count_transform}{Count normalization transform passed to
\code{\link{scscale_fit}}.}
\item{mp_max_iter, mp_grid_n}{MAD/MP spike-fitting controls.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{...}{Additional arguments passed to \code{\link{scscale_fit}}.}
}
\description{
Fit both modalities, store their normalized matrices, compute right singular
vectors, and align the low-rank subspaces.
}
\value{
A \code{scscale_pair_fit} object with \code{x_fit}, \code{y_fit},
\code{z_X}, \code{z_Y}, the overlap matrix \code{P}, baseline \code{mi}, and
\code{I_infinity}.
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts[1:100, 1:250]
adt <- gse164378_3p_citeseq_hvg$adt_counts[1:35, 1:250]
pair <- scscale_pair_fit(rna, adt, mp_grid_n = 600, mp_max_iter = 50)
pair$I_infinity
}
\keyword{models}
