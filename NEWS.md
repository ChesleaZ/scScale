# scScale 0.0.4

- `scscale_effective_fit()` now constrains `I_infinity` to at most twice the
  maximum empirical MI by default. Fits reaching this ceiling retain their
  observed-range curve but set `I_infinity_identified = FALSE`, record the
  upper boundary, and emit an explicit warning. Set
  `I_infinity_max_factor = Inf` to disable this safeguard.
- The UMI-anchored effective model now applies a parameter-free shifted
  softplus to the linear signal, keeping it strictly above `sqrt(d / n)`,
  before evaluating the original RMT recoverability. This removes the
  residual threshold transition while retaining the linear strong-signal
  limit. The UMI-anchored cell-number law holds
  those fitted UMI parameters fixed and estimates only its scale \(\gamma\),
  without applying a BBP indicator or an additional gate to the correction.
  Its fixed reference level is now the empirical MI at \((n_0,U_0)\), rather
  than the possibly biased theoretical UMI-fit value; this anchor is observed
  and therefore does not add a fitted parameter.
- `scscale_empirical_mi()` and `scscale_empirical_inu_grid()` now accept
  separate `x_rank` and `y_rank` values for matrix--matrix MI while retaining
  `r` as the backward-compatible shared-rank default.
- The canonical UMI calibration now constrains every component law to
  `q_X,k(U) = b_k U`. Enforcing `q_X,k(0) = 0` removes the large positive
  intercepts that previously kept recoverability nearly saturated at low
  UMI depth.
- Calibrated spike-model MI now deattenuates the observed reference PC
  alignment before extrapolation. The fitted object stores the observed
  overlap as `alignment$P` and the regularized calibration as `alignment$K`;
  predictions apply state-specific recoverabilities to `K`, avoiding the
  previous double attenuation at the reference dataset.
- Changed the canonical UMI spike calibration to estimate
  \(q_{X,k}(U)=[\widehat v_k(U)-\widehat\tau_X^2(U)]_+/
  \widehat\tau_X^2(U)\) directly along fixed reference directions. The noise
  estimate now comes from the fixed full-depth residual subspace, and the
  package fits \(q_{X,k}(U)=b_kU\) by zero-anchored ordinary least squares.
- Cell-number scaling is now a pure prediction: \(q_{X,k}(U)\) is independent
  of \(n\), cell number enters only through \(c_X(n)=p_X/n\), and the package
  then calculates \(\theta_X(n,U)\) and MI.
- Added `scscale_spike_fit()` as the canonical complete spike-model
  calibration, with stable reference directions, per-spike UMI parameters,
  and one shared state evaluator.
- Added fit-object interfaces for `scscale_cell_scaling()`,
  `scscale_umi_scaling()`, and `scscale_joint_scaling()`.
- Added `scscale()` as the single user-facing entry point for the complete
  paired UMI-by-cell-number scaling workflow.
- Added `plot.scscale_model()` for UMI, cell-number, and joint-grid scaling
  plots.
- UMI plots now use a dense curve reconstructed from the fitted per-spike
  linear UMI parameters instead of connecting raw depth refits.

# scScale 0.0.2

- `scscale_umi_mi()` now saves a long per-reference-spike `scaling` table.
- Added per-spike linear UMI parameters in `q_fit_rate`: intercept, slope,
  fitted R-squared, number of sampled rates, and number of detected rates.
- Added `q_X_rate_hat` to the long table and total-spike summaries in `q_total`
  and `q_total_rate_fit`.
# scScale 0.0.5

- Added a common representation API for external embeddings and the standard
  PCA preprocessing path.
- Added UMI-anchored and held-out cell-number fitting wrappers that share the
  same empirical MI and effective-model plumbing.
