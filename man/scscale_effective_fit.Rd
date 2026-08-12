\name{scscale_effective_fit}
\alias{scscale_effective_fit}
\alias{scscale_effective_cell_calibrate}
\alias{scscale_effective_cell_calibrate_B}
\alias{scscale_effective_mi}
\alias{scscale_effective_recoverability}
\alias{scscale_effective_signal}
\alias{predict.scscale_effective_fit}
\alias{print.scscale_effective_fit}
\title{Fit a UMI-Anchored Effective Scaling Model}
\description{
Fits the five-parameter effective model
\eqn{(I_\infty,C,\theta_Y,a,b)} to empirical UMI-scaling mutual information.
The nonlinear optimizer varies only \eqn{(\theta_Y,a,b)}. At every objective
evaluation, \eqn{(I_\infty,C)} are solved by constrained weighted least
squares. An optional sixth parameter, \eqn{\gamma}, calibrates cell-number
scaling around the fixed reference count without changing the UMI fit.
}
\usage{
scscale_effective_fit(
  data,
  d,
  n_reference,
  U_reference = NULL,
  U_col = "U",
  mi_col = "I_empirical",
  weighting = c("inverse_variance", "equal"),
  start = NULL,
  n_starts = 1L,
  seed = 1,
  I_infinity_max_factor = 2,
  control = list(maxit = 5000)
)

scscale_effective_recoverability(q, n, d)

scscale_effective_signal(q_linear, n, d)

scscale_effective_cell_calibrate(
  object,
  data,
  n_col = "n",
  U_col = "U",
  mi_col = "I_empirical",
  weighting = c("equal", "inverse_variance"),
  gamma_bounds = c(0, 1e4)
)

scscale_effective_cell_calibrate_B(
  object,
  data,
  n_col = "n",
  U_col = "U",
  mi_col = "I_empirical",
  weighting = c("equal", "inverse_variance"),
  B_bounds = c(1e-12, 1e12)
)

scscale_effective_mi(
  n, U, d, I_infinity, C, theta_Y, a, b
)

\method{predict}{scscale_effective_fit}(
  object,
  newdata = NULL,
  n = NULL,
  U = NULL,
  type = c("response", "all"),
  ...
)

\method{print}{scscale_effective_fit}(x, ...)
}
\arguments{
  \item{data}{A data frame containing UMI levels and empirical MI observations.
  For cell calibration, it contains cell-number and empirical MI observations.
  Replicate rows at the same level are summarized automatically.}
  \item{d}{Fixed representation dimension.}
  \item{n_reference}{Cell number used for the UMI calibration experiment.}
  \item{U_reference}{Reference depth used for cell-number prediction. Defaults
  to the largest observed UMI level.}
  \item{U_col,mi_col}{Names of the UMI and empirical MI columns.}
  \item{weighting}{Use inverse variance of replicate means when available, or
  equal weights.}
  \item{n_col}{Name of the cell-number column.}
  \item{gamma_bounds}{Two nonnegative values bounding the one-dimensional
  optimization of \eqn{\gamma}.}
  \item{B_bounds}{Two increasing positive values bounding the direct
  reference-scale optimization.}
  \item{start}{Optional list containing starting values for
  \code{theta_Y}, \code{a}, and \code{b}.}
  \item{n_starts}{Number of nonlinear optimization starts. One is the default.}
  \item{seed}{Random seed used only when \code{n_starts > 1}.}
  \item{I_infinity_max_factor}{Upper bound for \eqn{I_\infty}, expressed as
  a multiple of the maximum empirical MI. The default is 2; use \code{Inf}
  to disable the upper bound. A fit reaching the bound is marked as not
  identified and emits a warning.}
  \item{control}{Control list passed to \code{optim}.}
  \item{q}{Linear effective signal strengths before applying the smooth floor.}
  \item{q_linear}{Linear effective signal strengths.}
  \item{n}{Positive cell numbers. For prediction, values are crossed with
  \code{U}.}
  \item{U}{Nonnegative median UMI counts per cell.}
  \item{I_infinity,C,theta_Y,a,b}{Effective-model parameters.}
  \item{object,x}{A fitted \code{scscale_effective_fit} object.}
  \item{newdata}{Optional data frame containing paired \code{n} and \code{U}
  columns.}
  \item{type}{Return fitted MI values or all intermediate quantities.}
  \item{...}{Unused.}
}
\details{
The linear signal, smoothly floored signal, and recoverability are
\deqn{q_{lin}(U)=a+bU, \qquad q_0(n)=\sqrt{d/n},}
\deqn{
q(n,U)=q_0(n)+q_0(n)\log\left[1+\exp\left(
\frac{q_{lin}(U)-q_0(n)}{q_0(n)}\right)\right],
}
\deqn{
\theta_X(n,U)=\frac{q(n,U)^2-d/n}{q(n,U)[q(n,U)+1]}.
}
The shifted softplus is monotone, stays strictly above \eqn{\sqrt{d/n}}, and
approaches the linear signal for strong signals. Its width is fixed by the
same RMT scale, so it introduces no additional fitted parameter or phase
transition.
The fitted mutual information is
\deqn{
I(n,U)=I_\infty-C\log\left[
\frac{1-\theta_X(n,U)\theta_Y}{1-\theta_Y}
\right].
}
Only UMI-scaling observations at \code{n_reference} enter the initial fit.
\code{scscale_effective_cell_calibrate()} holds those five parameters fixed and
estimates only \eqn{\gamma} from cell-number observations, using
\deqn{
A_n(U)=\frac{\theta_Y d}
{q(U)\{q(U)+1\}\{1-\theta_Yq(U)/(q(U)+1)\}},
}
and replacing \eqn{A_n(U)} by \eqn{\gamma A_n(U)} in the cell-number
law. With \eqn{I_0} denoting the empirical MI at the observed reference point
\eqn{(n_0,U_0)}, the prediction is
\deqn{
I_\gamma(n,U)=I_0+C
\log\left\{\frac{1+\gamma A_n(U)/n_0}{1+\gamma A_n(U)/n}\right\}.
}
The five UMI-fit parameters are held fixed and only \eqn{\gamma} is estimated.
The empirical anchor is observed rather than fitted, so it adds no free
parameter. No BBP indicator or additional gate is applied to the cell-number
correction.
}
\value{
\code{scscale_effective_fit()} returns an object containing coefficients,
the summarized calibration data, fit metrics, and optimization diagnostics.
\code{scscale_effective_cell_calibrate()} returns the same object with
\code{gamma}, the empirical and theoretical values in \code{cell_anchor},
cell calibration data, and cell fit metrics added.
\code{scscale_effective_cell_calibrate_B()} instead estimates one nonnegative
direct reference scale, \code{B_reference}, while retaining the empirical
anchor and fixed UMI-fit parameters.
\code{predict()} returns fitted MI values or a data frame containing
\code{n}, \code{U}, \code{q_linear}, \code{q}, \code{theta_X}, \code{A_n}, \code{n_BBP},
and \code{I_fitted}.
}
\examples{
U <- c(100, 200, 400, 800, 1600, 3200)
truth <- scscale_effective_mi(
  n = 2000, U = U, d = 10,
  I_infinity = 1.5, C = 0.7, theta_Y = 0.6,
  a = 0.1, b = 0.001
)
fit <- scscale_effective_fit(
  data.frame(U = U, I_empirical = truth),
  d = 10,
  n_reference = 2000
)
predict(fit, n = c(500, 1000, 2000), U = max(U))
}
