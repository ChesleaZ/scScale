\name{scscale_pair_fit}
\alias{scscale_pair_fit}
\title{Fit Matched Modalities for MI Scaling}
\usage{
scscale_pair_fit(x, y, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), mp_max_iter = 300, mp_grid_n = 3000,
  y_mode = c("spike", "subspace", "auto"), y_rank = NULL,
  y_center = TRUE, y_scale = FALSE, y_input = c("auto", "matrix",
  "labels", "numeric"), auto_low_rank_max = 50, rank_tol = 1e-08,
  x_rank = NULL, x_force_spikes = FALSE, use_irlba = TRUE, ...)
}
\arguments{
\item{x}{Matched feature-by-cell count matrix for the noisy input modality.}
\item{y}{Matched target. Usually a feature-by-cell count matrix for a second
noisy modality, but can also be a deterministic low-rank target matrix or
label/numeric vector when \code{y_mode = "subspace"}.}
\item{count_transform}{Count normalization transform passed to
\code{\link{scscale_fit}}.}
\item{mp_max_iter, mp_grid_n}{MAD/MP spike-fitting controls.}
\item{y_mode}{How to handle \code{y}. \code{"spike"} preserves the original
paired-modality behavior and fits a spike model to \code{y}. \code{"subspace"}
treats \code{y} as a deterministic target subspace with target recoverability
equal to one. \code{"auto"} uses conservative low-rank/one-hot detection.}
\item{y_rank}{Optional target rank when \code{y_mode = "subspace"}. Defaults to
the numerical rank after preprocessing.}
\item{y_center, y_scale}{Whether to center or scale deterministic target rows
before taking right singular vectors.}
\item{y_input}{Target type for deterministic targets. \code{"auto"} infers
matrix, label, or numeric-vector input.}
\item{auto_low_rank_max, rank_tol}{Controls for conservative automatic
deterministic-target detection.}
\item{x_rank}{Optional number of input spikes to retain.}
\item{x_force_spikes}{Whether to force the retained input components to be
treated as spikes.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{...}{Additional arguments passed to \code{\link{scscale_fit}}.}
}
\description{
Fit the noisy input modality, optionally fit or directly construct the target
side, compute right singular vectors, and align the low-rank subspaces.
}
\value{
A \code{scscale_pair_fit} object with \code{x_fit}, \code{y_fit},
\code{z_X}, \code{z_Y}, the overlap matrix \code{P}, baseline \code{mi}, and
\code{I_infinity}. When \code{y_mode = "subspace"}, \code{y_fit} is a
\code{scscale_subspace_fit} object and no MP spike model is fitted to the
target.
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts[1:100, 1:250]
adt <- gse164378_3p_citeseq_hvg$adt_counts[1:35, 1:250]
pair <- scscale_pair_fit(rna, adt, mp_grid_n = 600, mp_max_iter = 50)
pair$I_infinity

labels <- factor(rep(c("A", "B"), each = 125))
label_pair <- scscale_pair_fit(
  rna,
  labels,
  y_mode = "subspace",
  y_input = "labels",
  mp_grid_n = 600,
  mp_max_iter = 50
)
label_pair$y_mode
}
\keyword{models}
