\name{scscale_empirical_fit_helpers}
\alias{scscale_empirical_values}
\alias{scscale_empirical_parameters}
\alias{scscale_predict_empirical_cell}
\alias{scscale_predict_empirical_umi}
\alias{scscale_predict_empirical_joint}
\alias{scscale_plot_empirical_fit}
\alias{scscale_plot_empirical_umi_levels}
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
}
\arguments{
\item{object}{A \code{scscale_empirical_inu_fit} object returned by
\code{\link{scscale_empirical_scaling_fit}}.}
\item{summarize}{Whether to aggregate empirical MI values over replicate
samples at each grid point.}
\item{cell_sampling_rate, sampling_rate}{Cell-number and UMI sampling rates.}
\item{n, U}{Absolute cell numbers and UMI depths. Provide rates or absolute
values, not both.}
\item{n_points, cell_points, umi_points}{Number of points used when generating
regular prediction curves.}
\item{type}{Which marginal fit to plot.}
\item{...}{Additional arguments passed to base \code{\link[graphics]{plot}}.}
}
\description{
Convenience helpers for working with a fitted direct empirical scaling object.
They expose the sampled empirical MI values, fitted parameters, marginal
predictions, joint predictions over the full cell-by-UMI combination, and base R
diagnostic plots.
}
\value{
\code{scscale_empirical_values()} returns a data frame of empirical grid values
and fitted joint predictions. \code{scscale_empirical_parameters()} returns a
data frame of fitted parameters and ceiling diagnostics. Prediction helpers
return data frames with inputs and \code{I_fit}. Plot helpers return
\code{object} invisibly.
}
\keyword{models}
