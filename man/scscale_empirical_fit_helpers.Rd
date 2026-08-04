\name{scscale_empirical_fit_helpers}
\alias{scscale_empirical_values}
\alias{scscale_empirical_parameters}
\alias{scscale_predict_empirical_cell}
\alias{scscale_predict_empirical_umi}
\alias{scscale_predict_empirical_joint}
\alias{scscale_plot_empirical_fit}
\alias{scscale_plot_empirical_surface}
\alias{scscale_plot_empirical_umi_levels}
\alias{scscale_empirical_inu_grid_corrected}
\title{Inspect, Predict, and Plot Direct Empirical Scaling Fits}
\usage{
scscale_empirical_values(object, summarize = FALSE)

scscale_empirical_parameters(object)

scscale_predict_empirical_cell(object, cell_sampling_rate = NULL, n = NULL,
  n_points = NULL)

scscale_predict_empirical_umi(object, sampling_rate = NULL, U = NULL,
  n_points = NULL)

scscale_predict_empirical_joint(object, cell_sampling_rate = NULL, n = NULL,
  sampling_rate = NULL, U = NULL, cell_points = NULL, umi_points = NULL)

scscale_plot_empirical_fit(object, type = c("both", "cell", "umi"),
  n_points = 200, ...)

scscale_plot_empirical_umi_levels(object, sampling_rate = NULL, U = NULL,
  n_points = 300, ...)

scscale_plot_empirical_surface(object, n = NULL, U = NULL, n_points = 60,
  U_points = 60, show_points = TRUE, theta = 35, phi = 25,
  expand = 0.65, surface_col = "#56B4E9", point_col = "#D55E00",
  border = NA, shade = 0.35, ...)

scscale_empirical_inu_grid_corrected(x, target, batch = NULL, n_grid = NULL,
  cell_sampling_rates = c(0.125, 0.1875, 0.25, 0.375, 0.50, 0.75,
  1.00), U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85,
  1.00), r = 10, cell_replicates = 10, umi_replicates = 1,
  n_workers = 1, correction = c("custom", "scvi"), correction_fun = NULL,
  correction_cache = NULL, correction_options = list(), seed = 1,
  input = c("counts", "normalized"), corrected_input = c("normalized",
  "counts"), target_input = input, target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE, scale = FALSE, eps = 1e-12, use_irlba = TRUE)
}
\arguments{
\item{object}{A \code{scscale_empirical_inu_fit} object returned by
\code{\link{scscale_empirical_scaling_fit}}.}
\item{x}{Feature-by-cell source matrix.}
\item{target}{Feature-by-cell target matrix or cell-level target vector.}
\item{batch}{Optional cell-level batch vector passed to \code{correction_fun}.}
\item{n_grid, cell_sampling_rates, U_grid, sampling_rates}{Cell-number and UMI
sampling grids. Provide one cell grid and one UMI grid.}
\item{r}{Fixed rank used for empirical subspace MI.}
\item{cell_replicates, umi_replicates}{Cell-sampling and full-matrix UMI
thinning replicates.}
\item{n_workers}{Number of parallel workers for cell-sampling replicates.}
\item{correction, correction_fun, correction_cache, correction_options}{Staged
correction controls. \code{correction_fun} receives a full UMI-thinned count
matrix, the aligned target, and the aligned batch vector, and must return a
feature-by-cell corrected matrix.}
\item{seed}{Random seed base.}
\item{input, corrected_input, target_input, target_depth, count_transform,
center, scale, eps, use_irlba}{Preprocessing and empirical MI controls.}
\item{summarize}{Whether to aggregate empirical MI values over replicate
samples at each grid point.}
\item{cell_sampling_rate, sampling_rate}{Cell-number and UMI sampling rates.}
\item{n, U}{Absolute cell numbers and UMI depths. Provide rates or absolute
values, not both.}
\item{n_points, U_points, cell_points, umi_points}{Number of points used when
generating regular prediction curves or the joint surface.}
\item{type}{Which marginal fit to plot.}
\item{show_points}{Whether to overlay the empirical grid means on the fitted
surface.}
\item{theta, phi, expand, surface_col, point_col, border, shade}{Graphical
controls passed to or used with \code{\link[graphics]{persp}}.}
\item{...}{Additional arguments passed to base \code{\link[graphics]{plot}}.}
}
\description{
Convenience helpers for working with a fitted direct empirical scaling object.
They expose the sampled empirical MI values, fitted parameters, marginal
predictions, joint predictions over the full cell-by-UMI combination, and base R
diagnostic plots, including a perspective surface of fitted mutual information
over absolute cell number and UMI depth.
\code{scscale_empirical_inu_grid_corrected()} is the staged
batch-aware grid engine used by \code{scscale_empirical_scaling_fit()} when
correction is requested: it UMI-thins the full raw matrix, applies a correction
function once per UMI level/replicate, and then computes MI over replicated cell
subsamples from the corrected representation.
}
\value{
\code{scscale_empirical_values()} returns a data frame of empirical grid values
and fitted joint predictions. \code{scscale_empirical_parameters()} returns a
data frame of fitted parameters and ceiling diagnostics. Prediction helpers
return data frames with inputs and \code{I_fit}. The marginal plot helpers
return \code{object} invisibly. \code{scscale_plot_empirical_surface()} returns
the prediction grid, perspective transformation matrix, and displayed
observations invisibly. \code{scscale_empirical_inu_grid_corrected()} returns
an empirical MI grid data frame.
}
\keyword{models}
