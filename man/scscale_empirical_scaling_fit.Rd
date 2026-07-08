\name{scscale_empirical_scaling_fit}
\alias{scscale_empirical_scaling_fit}
\alias{scscale_fit_empirical_cell_law}
\alias{scscale_fit_empirical_umi_law}
\alias{scscale_fit_empirical_joint_law}
\alias{scscale_empirical_global_ceiling}
\alias{scscale_empirical_cell_curve}
\alias{scscale_empirical_umi_curve}
\alias{scscale_empirical_joint_curve}
\title{Fit Direct Empirical Scaling Laws}
\usage{
scscale_empirical_scaling_fit(x, target, n_grid, U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85,
  1.00), r = 10, n_replicates = 10, n_workers = 1, n_ref = NULL,
  U_ref = NULL, joint_mode = c("hybrid", "transfer", "free"),
  seed = 1, input = c("counts", "normalized"), target_input = input,
  target_depth = 1e4, count_transform = c("log1p_cpm",
  "pearson_residual", "log1p"), center = TRUE, scale = FALSE,
  eps = 1e-12, use_irlba = TRUE, cell_start = NULL, umi_start = NULL,
  joint_start = NULL, cell_control = list(maxit = 10000),
  umi_control = list(maxit = 10000),
  joint_control = list(maxit = 20000))

scscale_fit_empirical_cell_law(data, n_col = "n", mi_col = "I_empirical",
  start = NULL, control = list(maxit = 10000))

scscale_fit_empirical_umi_law(data, U_col = "U", mi_col = "I_empirical",
  start = NULL, control = list(maxit = 10000))

scscale_fit_empirical_joint_law(data, n_col = "n", U_col = "U",
  mi_col = "I_empirical", cell_fit = NULL, umi_fit = NULL,
  mode = c("hybrid", "transfer", "free"), start = NULL,
  control = list(maxit = 20000))

scscale_empirical_global_ceiling(cell_fit, umi_fit, n0, U0)

scscale_empirical_cell_curve(n, I_infinity, C_n, A_n)

scscale_empirical_umi_curve(U, I_infinity, C_U, A_U, U_offset)

scscale_empirical_joint_curve(n, U, I_infinity, C_n, A_n, C_U, A_U,
  U_offset)
}
\arguments{
\item{x}{Feature-by-cell source matrix.}
\item{target}{Feature-by-cell target matrix or cell-level target vector.}
\item{n_grid}{Cell numbers to sample.}
\item{U_grid}{Target UMI depths per cell. Provide either \code{U_grid} or
\code{sampling_rates}. If \code{U_grid} is supplied and \code{sampling_rates}
is omitted, the default sampling-rate grid is ignored.}
\item{sampling_rates}{Count downsampling fractions in \code{(0,1]}.}
\item{r}{Fixed rank used for empirical subspace MI.}
\item{n_replicates}{Number of replicate subsamples per grid point.}
\item{n_workers}{Number of parallel workers for empirical grid evaluation.
The default \code{1} evaluates sequentially.}
\item{n_ref, U_ref}{Reference cell number and UMI depth for marginal fits.
When omitted, the largest sampled value is used. If a supplied reference is not
on the grid, the nearest sampled grid value is used.}
\item{joint_mode}{Joint fitting mode passed to
\code{scscale_fit_empirical_joint_law}.}
\item{seed}{Random seed base for grid sampling.}
\item{input, target_input}{Input type for \code{x} and matrix targets.}
\item{target_depth, count_transform, center, scale}{Preprocessing controls
passed to empirical MI calculation.}
\item{eps}{Numerical clipping value for MI.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial SVD when available.}
\item{cell_start, umi_start, joint_start}{Optional named starting values for
the marginal and joint nonlinear fits.}
\item{cell_control, umi_control, joint_control}{Control lists passed to
\code{\link[stats]{optim}} for the marginal and joint nonlinear fits.}
\item{data}{Data frame containing empirical MI observations.}
\item{n_col, U_col, mi_col}{Column names for cell number, UMI depth, and empirical MI.}
\item{start}{Optional named list of starting values.}
\item{control}{Control list passed to \code{\link[stats]{optim}}.}
\item{cell_fit, umi_fit}{Marginal empirical scaling fits used to initialize or transfer joint-law parameters.}
\item{mode}{Joint fitting mode. \code{"transfer"} fixes marginal gap parameters and fits only the global ceiling; \code{"hybrid"} fixes scale parameters and refits the ceiling and amplitudes; \code{"free"} refits all joint parameters.}
\item{n, U}{Cell number and UMI depth values for curve evaluation.}
\item{n0, U0}{Reference cell number and UMI depth for extrapolating the global ceiling from marginal fits.}
\item{I_infinity, C_n, A_n, C_U, A_U, U_offset}{Scaling-law parameters.}
}
\description{
Fit direct nonlinear scaling laws to fixed-rank empirical subspace MI values,
without spike calling or Marchenko--Pastur bulk fitting. The high-level
\code{scscale_empirical_scaling_fit()} function samples an empirical
\code{I(n,U)} grid, preprocesses each sampled view, computes empirical MI, and
fits the marginal and joint laws into one \code{scscale_empirical_inu_fit}
object. The cell law is
\deqn{I(n)=I_{\infty|U_0}-C_n\log(1+A_n/n).}
The UMI law is
\deqn{I(U)=I_{\infty|n_0}-C_U\log(1+A_U/(U+U_0)).}
The joint law is additive in the two empirical information gaps:
\deqn{I(n,U)=I_\infty-C_n\log(1+A_n/n)-C_U\log(1+A_U/(U+U_0)).}
}
\value{
\code{scscale_empirical_scaling_fit()} returns an object of class
\code{scscale_empirical_inu_fit}, containing the empirical grid, marginal fits,
joint fit, reference values, ceiling diagnostics, residuals, and fit metrics.
The lower-level fitting functions return objects of class
\code{scscale_empirical_scaling_fit}, with fitted coefficients, residuals, fit
metrics, and \code{predict()} support. \code{scscale_empirical_global_ceiling}
returns the two marginal extrapolations to \code{I(infinity, infinity)} and
their average.
}
\keyword{models}
