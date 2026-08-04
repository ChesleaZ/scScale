\name{scscale_spike_workflow}
\alias{scscale_spike_fit}
\alias{scscale_cell_scaling}
\alias{scscale_umi_scaling}
\alias{scscale_joint_scaling}
\alias{print.scscale_spike_fit}
\title{Calibrate and Evaluate the Complete Gaussian Spike Model}
\usage{
scscale_spike_fit(x, y,
  umi_sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  umi_replicates = 1,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  mp_max_iter = 300, mp_grid_n = 3000, x_rank = NULL, y_rank = NULL,
  x_force_spikes = FALSE, y_force_spikes = FALSE,
  alignment_theta_floor = 1e-6,
  alignment_max_singular = 1 - 1e-8, seed = 1,
  use_irlba = TRUE, ...)

scscale_cell_scaling(fit = NULL, n = NULL, ..., q_X = NULL)

scscale_umi_scaling(fit = NULL, sampling_rate = NULL, U = NULL, ...,
  counts = NULL)

scscale_joint_scaling(fit, n, sampling_rate = NULL, U = NULL,
  eps = 1e-12)

\method{print}{scscale_spike_fit}(x, ...)
}
\arguments{
\item{x, y}{Matched feature-by-cell count matrices for the X and Y
modalities.}
\item{umi_sampling_rates}{UMI thinning fractions used to calibrate the
per-component linear \eqn{q_X(U)} model. Each fraction is converted to an
absolute depth using the reference median UMI count per cell. The full-depth
rate 1 is added automatically when absent.}
\item{umi_replicates}{Number of independent count-thinning replicates per UMI
rate.}
\item{count_transform}{Count normalization transform.}
\item{mp_max_iter, mp_grid_n}{MP bulk and spike-selection controls.}
\item{x_rank, y_rank}{Optional maximum numbers of retained X and Y spike
components.}
\item{x_force_spikes, y_force_spikes}{Whether to force the retained leading
components to be spikes.}
\item{alignment_theta_floor}{Lower bound applied to reference recoverabilities
when deattenuating the observed reference alignment matrix.}
\item{alignment_max_singular}{Upper bound for singular values of the
deattenuated alignment matrix.}
\item{seed}{Random seed base for UMI thinning.}
\item{use_irlba}{Whether to use \pkg{irlba} for partial singular-value
decompositions when available.}
\item{...}{Additional preprocessing and spike-fit arguments. The scaling
functions retain compatibility with their earlier low-level signatures.}
\item{fit}{A complete \code{scscale_spike_fit} object.}
\item{n}{One or more absolute cell numbers.}
\item{sampling_rate}{One or more positive UMI sampling rates. Values above one
request extrapolation.}
\item{U}{Optional absolute median UMI depths per cell, converted to sampling
rates using the reference median depth. Supply \code{sampling_rate} or
\code{U}, not both.}
\item{eps}{Numerical clipping value for the MI calculation.}
\item{q_X, counts}{Compatibility arguments for the lower-level cell and
UMI helpers.}
}
\description{
\code{scscale_spike_fit()} performs all data-dependent calibration once. It
fits the reference X and Y Gaussian spike models, selects ranks, computes
reference \eqn{q} and recoverability parameters, measures the observed X--Y
PC alignment \eqn{P_{\rm obs}}, projects every UMI-thinned X matrix onto the fixed full-depth feature
directions, estimates the noise from their fixed residual complement, and
calculates
\deqn{q_{X,k}(U)=
\frac{[\widehat v_k(U)-\widehat\tau_X^2(U)]_+}
{\widehat\tau_X^2(U)}.}
It then fits a zero-anchored ordinary least-squares UMI law for every reference
X spike:
\deqn{q_{X,k}(U)=b_kU.}
The zero-depth constraint enforces \eqn{q_{X,k}(0)=0} and prevents a fitted
intercept from creating nonzero recoverability without any sampled UMIs.

The three scaling helpers are deterministic views of the same fitted model.
\code{scscale_cell_scaling()} evaluates \eqn{I(n,1)},
\code{scscale_umi_scaling()} evaluates \eqn{I(n_{\mathrm{ref}},\rho)}, and
\code{scscale_joint_scaling()} evaluates the Cartesian grid
\eqn{I(n,\rho)}. Cell number does not enter \eqn{q_{X,k}(U)}. Each helper
calculates \eqn{c_X(n)=p_X/n}, applies the RMT recoverability formula to
\eqn{q_{X,k}(U)} and \eqn{c_X(n)}, and then computes MI using the fixed
target-side recoverability and a calibrated alignment matrix. The calibration
removes reference recoverability from the observed PC overlap,
\deqn{\widehat K=\Theta_{X,{\rm ref}}^{-1/2}P_{\rm obs}
\Theta_{Y,{\rm ref}}^{-1/2},}
then clips singular values of \eqn{\widehat K} to the configured upper bound.
At a requested state it reconstructs
\deqn{M(n,U)=\Theta_X(n,U)^{1/2}\widehat K
\Theta_Y(n,U)^{1/2}.}
When clipping is inactive, this reproduces the empirical reference overlap
exactly instead of attenuating it twice.

The spike-model ceiling is calculated once from the fitted object by setting
every retained X-side recoverability to one while leaving the fitted
target-side recoverabilities unchanged:
\deqn{I_\infty=
-\frac{1}{2}\sum_\ell
\log\left[1-\sigma_\ell^2\left(\widehat K\Theta_Y^{1/2}\right)\right].}
It is therefore a derived quantity rather than a fitted parameter, and the
same scalar is returned by the cell-number, UMI-depth, and joint scaling
helpers.
}
\value{
\code{scscale_spike_fit()} returns an object of class
\code{scscale_spike_fit}. Its \code{x} and \code{y} elements contain the
canonical component parameter tables; \code{alignment$P} contains the
observed reference overlap matrix and \code{alignment$K} contains its
deattenuated, singular-value-regularized calibration; \code{umi$parameters}
contains the zero intercepts and fitted \eqn{q} slopes; and \code{diagnostics} contains the
rate-specific projected variance, excess variance, fixed-residual noise
estimates, and reference MP histories.

The scaling helpers return data frames containing the requested inputs,
aspect ratios, theoretical MI, the shared \eqn{I_\infty} upper bound, and
effective rank. The underlying per-component states are attached to these
data frames as a \code{states} attribute. The fitted object also stores the
scalar as \code{I_infinity} and its construction in \code{ceiling}.
}
\examples{
\dontrun{
data(gse164378_3p_citeseq_hvg)
fit <- scscale_spike_fit(
  gse164378_3p_citeseq_hvg$rna_counts,
  gse164378_3p_citeseq_hvg$adt_counts
)
cell <- scscale_cell_scaling(fit, n = c(500, 1000, 2000, 4000))
umi <- scscale_umi_scaling(fit, sampling_rate = seq(0.1, 1, length.out = 50))
joint <- scscale_joint_scaling(
  fit,
  n = c(500, 1000, 2000, 4000),
  sampling_rate = c(0.2, 0.5, 1)
)
}
}
\seealso{
\code{\link{scscale_fit}},
\code{\link{scscale_empirical_scaling_fit}}
}
\keyword{models}
