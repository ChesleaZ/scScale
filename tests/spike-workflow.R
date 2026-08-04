library(scScale)

project_reference_spikes <- getFromNamespace(
  "scscale_project_reference_spikes",
  "scScale"
)
fixed_residual_tau2 <- getFromNamespace(
  "scscale_fixed_residual_tau2",
  "scScale"
)
fit_umi_q_parameters <- getFromNamespace(
  "scscale_fit_umi_q_parameters",
  "scScale"
)
spike_state <- getFromNamespace("scscale_spike_state", "scScale")
calibrate_reference_alignment <- getFromNamespace(
  "scscale_calibrate_reference_alignment",
  "scScale"
)

matrix_uncentered <- matrix(c(3, 1, 4, 1, 5, 9), nrow = 2)
attr(matrix_uncentered, "row_center") <- rowMeans(matrix_uncentered)
reference_vectors <- diag(2)
projection <- project_reference_spikes(
  matrix_uncentered,
  reference_vectors,
  tau2 = 2
)
matrix_centered <- sweep(
  matrix_uncentered,
  1,
  rowMeans(matrix_uncentered),
  "-"
)
expected_variance <- rowSums(matrix_centered^2) / ncol(matrix_centered)
stopifnot(isTRUE(all.equal(
  projection$projected_variance,
  expected_variance
)))
stopifnot(isTRUE(all.equal(
  projection$q_X,
  pmax(expected_variance - 2, 0) / 2
)))

residual_matrix <- rbind(
  c(2, -2, 2, -2),
  c(1, -1, -1, 1),
  c(3, -3, -3, 3)
)
stopifnot(isTRUE(all.equal(
  fixed_residual_tau2(residual_matrix, diag(3)[, 1, drop = FALSE]),
  sum(residual_matrix[-1, , drop = FALSE]^2) / (4 * 2)
)))

calibration <- data.frame(
  component = rep(1:2, each = 3),
  U = rep(c(10, 20, 30), 2),
  q_X = c(3, 5, 7, 4, 7, 10)
)
parameters <- suppressWarnings(
  fit_umi_q_parameters(calibration, reference_U = 30)
)
stopifnot(
  isTRUE(all.equal(parameters$q_intercept_a, c(0, 0))),
  isTRUE(all.equal(
    parameters$q_slope_b,
    c(17 / 70, 12 / 35)
  )),
  isTRUE(all.equal(
    parameters$q_fit_reference,
    30 * c(17 / 70, 12 / 35)
  ))
)

P_observed <- diag(c(0.4, 0.2))
theta_X_reference <- c(0.8, 0.5)
theta_Y_reference <- c(0.7, 0.6)
alignment_calibration <- calibrate_reference_alignment(
  P_observed,
  theta_X = theta_X_reference,
  theta_Y = theta_Y_reference
)
P_reconstructed <- sweep(
  alignment_calibration$K,
  1L,
  sqrt(theta_X_reference),
  "*"
)
P_reconstructed <- sweep(
  P_reconstructed,
  2L,
  sqrt(theta_Y_reference),
  "*"
)
stopifnot(
  isTRUE(all.equal(P_reconstructed, P_observed, tolerance = 1e-12)),
  alignment_calibration$reference_reconstruction_error < 1e-12,
  !any(alignment_calibration$singular_clipped)
)

fit <- structure(
  list(
    reference = list(
      n = 100,
      p_X = 20,
      p_Y = 10,
      median_umi_per_cell_X = 30
    ),
    x = list(
      parameters = data.frame(
        q_intercept_a = 1,
        q_slope_b = 0.2
      )
    ),
    y = list(
      parameters = data.frame(
        q_ref = 2,
        theta_ref = 0.7,
        theta_infinity_ref = 0.8
      )
    ),
    alignment = list(
      P = matrix(0.5, nrow = 1, ncol = 1),
      K = matrix(0.6, nrow = 1, ncol = 1)
    )
  ),
  class = "scscale_spike_fit"
)
state_50 <- spike_state(fit, n = 50, sampling_rate = 0.5)
state_100 <- spike_state(fit, n = 100, sampling_rate = 0.5)
state_full_depth <- spike_state(fit, n = 100, sampling_rate = 1)
expected_state_50_mi <- scscale_low_rank_mi(
  theta_X = state_50$theta_X,
  theta_Y = 0.7,
  P = matrix(0.6, nrow = 1, ncol = 1)
)$mi
expected_I_infinity <- scscale_low_rank_mi(
  theta_X = 1,
  theta_Y = 0.7,
  P = matrix(0.6, nrow = 1, ncol = 1)
)$mi
stopifnot(
  isTRUE(all.equal(state_50$U, 15)),
  isTRUE(all.equal(state_50$q_X, 3)),
  isTRUE(all.equal(state_100$q_X, state_50$q_X)),
  state_100$theta_X > state_50$theta_X,
  state_100$I_theory > state_50$I_theory,
  isTRUE(all.equal(state_50$I_theory, expected_state_50_mi)),
  isTRUE(all.equal(state_50$theta_Y, 0.7)),
  isTRUE(all.equal(state_50$theta_X_infinity, 1)),
  isTRUE(all.equal(state_50$theta_Y_infinity, 0.7)),
  isTRUE(all.equal(state_50$I_infinity, expected_I_infinity)),
  isTRUE(all.equal(state_100$I_infinity, expected_I_infinity)),
  isTRUE(all.equal(state_full_depth$I_infinity, expected_I_infinity)),
  is.null(state_50$d2_X)
)
