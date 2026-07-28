\name{scscale}
\alias{scscale}
\alias{plot.scscale_model}
\title{Fit the Complete scScale Joint Model}
\usage{
scscale(
  x,
  y,
  umi_sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  n_grid = NULL,
  n_min = NULL,
  n_points = 12L,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  mp_max_iter = 300,
  mp_grid_n = 3000,
  seed = 1,
  use_irlba = TRUE,
  ...
)
\method{plot}{scscale_model}(
  x,
  type = c("umi", "cell", "joint"),
  ...
)
}
\arguments{
\item{x}{Feature-by-cell count matrix for the noisy input modality.}
\item{y}{Matched target matrix, label vector, or numeric target.}
\item{umi_sampling_rates}{Global UMI downsampling fractions in \code{(0,1]}.
These are the depths at which the spectral model is refitted.}
\item{n_grid}{Optional absolute cell-number grid. When omitted, a log-spaced
grid is generated automatically.}
\item{n_min}{Optional minimum absolute cell number for the automatic grid.
The default is \code{max(4 * number_of_X_spikes, 500)}, with a small-dataset
fallback.}
\item{n_points}{Number of points in the automatically generated cell grid.}
\item{count_transform}{Count normalization transform.}
\item{mp_max_iter, mp_grid_n}{MP bulk and spike-fitting controls.}
\item{seed}{Random seed for UMI downsampling.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{...}{For \code{scscale()}, additional arguments passed to
\code{\link{scscale_pair_fit}}, including target-mode and rank controls. For
the plot method, additional graphical arguments passed to \code{\link{plot}}.}
\item{type}{Scaling result to plot: UMI depth, cell number, or the joint
cell-number by UMI grid.}
}
\description{
Fit the full-depth spectral relationship between matched \code{x} and
\code{y}, learn UMI-depth scaling, calculate cell-number scaling, and return
the joint theoretical mutual-information grid \eqn{I(n,U)} in one object.
}
\value{
An \code{scscale_model} object containing \code{pair}, \code{umi},
\code{cell}, and \code{joint} results, along with the evaluated grids and
consolidated parameter tables.
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts[1:100, 1:250]
adt <- gse164378_3p_citeseq_hvg$adt_counts[1:35, 1:250]
model <- scscale(
  rna,
  adt,
  umi_sampling_rates = c(0.5, 1),
  n_grid = c(125, 250),
  mp_grid_n = 600,
  mp_max_iter = 50
)
model
plot(model, type = "umi")
plot(model, type = "cell")
plot(model, type = "joint")
}
\seealso{
\code{\link{scscale_pair_fit}},
\code{\link{scscale_umi_mi}},
\code{\link{scscale_cell_number_mi}},
\code{\link{scscale_cell_number_by_umi_mi}}
}
\keyword{models}
