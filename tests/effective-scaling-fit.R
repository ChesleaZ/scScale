library(scScale)

U <- c(50, 100, 200, 400, 800, 1600, 3200)
n_reference <- 2000
d <- 10
truth <- c(
  I_infinity = 1.4,
  C = 0.6,
  theta_Y = 0.65,
  a = 0.08,
  b = 0.0012
)

theta <- scscale_effective_recoverability(
  q = truth[["a"]] + truth[["b"]] * U,
  n = n_reference,
  d = d
)
q_threshold <- sqrt(d / n_reference)
q_linear_check <- q_threshold * c(-10, -2, -1, 0, 1, 2, 20)
q_gate_check <- scscale_effective_signal(
  q_linear = q_linear_check,
  n = n_reference,
  d = d
)
theta_gate_check <- scscale_effective_recoverability(
  q = q_linear_check,
  n = n_reference,
  d = d
)
theta_gate_expected <- (q_gate_check^2 - d / n_reference) /
  (q_gate_check * (q_gate_check + 1))
stopifnot(
  all(theta >= 0),
  all(theta <= 1),
  all(diff(theta) >= 0),
  all(q_gate_check > q_threshold),
  all(diff(q_gate_check) > 0),
  q_gate_check[[1]] - q_threshold < 1e-5,
  abs(q_gate_check[[7]] - q_linear_check[[7]]) < 1e-8,
  all(theta_gate_check > 0),
  all(diff(theta_gate_check) >= 0),
  isTRUE(all.equal(theta_gate_check, theta_gate_expected, tolerance = 1e-14))
)

mi <- scscale_effective_mi(
  n = n_reference,
  U = U,
  d = d,
  I_infinity = truth[["I_infinity"]],
  C = truth[["C"]],
  theta_Y = truth[["theta_Y"]],
  a = truth[["a"]],
  b = truth[["b"]]
)
stopifnot(all(diff(mi) >= 0), max(mi) <= truth[["I_infinity"]])

calibration <- data.frame(
  U = rep(U, each = 2),
  I_empirical = rep(mi, each = 2) + rep(c(-1, 1) * 1e-4, length(U))
)
fit <- scscale_effective_fit(
  calibration,
  d = d,
  n_reference = n_reference,
  U_reference = max(U),
  weighting = "inverse_variance",
  n_starts = 5,
  seed = 11
)
stopifnot(
  inherits(fit, "scscale_effective_fit"),
  identical(names(fit$coefficients), c("I_infinity", "C", "theta_Y", "a", "b")),
  fit$coefficients[["I_infinity"]] >= max(fit$calibration$I_mean),
  fit$coefficients[["C"]] >= 0,
  fit$coefficients[["theta_Y"]] > 0,
  fit$coefficients[["theta_Y"]] < 1,
  fit$coefficients[["a"]] >= 0,
  fit$coefficients[["b"]] > 0,
  fit$metrics$rmse < 0.01
)

linear_data <- data.frame(
  U = rep(U, each = 2),
  I_empirical = rep(seq(0.1, 0.7, length.out = length(U)), each = 2)
)
linear_fit <- suppressWarnings(scscale_effective_fit(
  linear_data,
  d = d,
  n_reference = n_reference,
  n_starts = 5,
  seed = 12,
  I_infinity_max_factor = 2
))
stopifnot(
  linear_fit$coefficients[["I_infinity"]] <=
    2 * max(linear_data$I_empirical) + 1e-8,
  is.logical(linear_fit$I_infinity_identified),
  linear_fit$I_infinity_boundary %in% c("none", "lower", "upper")
)

umi_prediction <- predict(fit, n = n_reference, U = U, type = "all")
cell_prediction <- predict(
  fit,
  n = c(250, 500, 1000, 2000, 4000),
  U = max(U),
  type = "all"
)
direct_prediction <- scscale_effective_mi(
  n = n_reference,
  U = U,
  d = d,
  I_infinity = fit$coefficients[["I_infinity"]],
  C = fit$coefficients[["C"]],
  theta_Y = fit$coefficients[["theta_Y"]],
  a = fit$coefficients[["a"]],
  b = fit$coefficients[["b"]]
)
stopifnot(
  nrow(umi_prediction) == length(U),
  nrow(cell_prediction) == 5L,
  all(diff(umi_prediction$I_fitted) >= -1e-10),
  all(diff(cell_prediction$I_fitted) >= -1e-10),
  isTRUE(all.equal(umi_prediction$I_fitted, direct_prediction, tolerance = 1e-10)),
  isTRUE(all.equal(
    predict(fit, newdata = cell_prediction[, c("n", "U")]),
    cell_prediction$I_fitted
  ))
)

gamma_truth <- 1.8
cell_n <- c(250, 500, 1000, 2000, 4000)
gamma_truth_fit <- fit
gamma_truth_fit$coefficients <- c(fit$coefficients, gamma = gamma_truth)
cell_mi <- predict(gamma_truth_fit, n = cell_n, U = max(U))
cell_data <- data.frame(
  n = rep(cell_n, each = 2),
  U = max(U),
  I_empirical = rep(cell_mi, each = 2) + rep(c(-1, 1) * 1e-5, length(cell_n))
)
calibrated <- scscale_effective_cell_calibrate(
  fit,
  cell_data,
  weighting = "equal"
)
uncalibrated_cell_rmse <- sqrt(mean(
  (cell_mi - predict(fit, n = cell_n, U = max(U)))^2
))
stopifnot(
  inherits(calibrated, "scscale_effective_fit"),
  "gamma" %in% names(calibrated$coefficients),
  !is.null(calibrated$cell_anchor),
  calibrated$coefficients[["gamma"]] >= 0,
  calibrated$cell_metrics$rmse < uncalibrated_cell_rmse,
  isTRUE(all.equal(
    predict(calibrated, n = n_reference, U = max(U)),
    calibrated$cell_anchor$I_empirical,
    tolerance = 1e-10
  ))
)

# The one-gamma cell law must be continuous through the diagnostic BBP point.
threshold <- calibrated$d /
  (calibrated$coefficients[["a"]] + calibrated$coefficients[["b"]] * max(U))^2
near_threshold <- predict(
  calibrated,
  n = threshold * c(1 - 1e-6, 1, 1 + 1e-6),
  U = max(U),
  type = "all"
)
stopifnot(
  max(abs(diff(near_threshold$I_fitted))) < 1e-5
)

# Direct-B calibration preserves both reference marginals and recovers a
# synthetic cell-number scale without depending on the absolute A(q) scale.
B_truth <- 750
cell_anchor_truth <- predict(fit, n = n_reference, U = max(U))
B_cell_mi <- cell_anchor_truth + fit$coefficients[["C"]] * (
  log1p(B_truth / n_reference) - log1p(B_truth / cell_n)
)
B_cell_data <- data.frame(
  n = rep(cell_n, each = 2),
  U = max(U),
  I_empirical = rep(B_cell_mi, each = 2) + rep(c(-1, 1) * 1e-5, length(cell_n))
)
B_calibrated <- scscale_effective_cell_calibrate_B(
  fit,
  B_cell_data,
  weighting = "equal"
)
stopifnot(
  inherits(B_calibrated, "scscale_effective_fit"),
  identical(B_calibrated$cell_mode, "B_reference"),
  isTRUE(all.equal(
    unname(B_calibrated$coefficients[["B_reference"]]),
    B_truth,
    tolerance = 1e-3
  )),
  B_calibrated$cell_metrics$rmse < 1e-4,
  isTRUE(all.equal(
    predict(B_calibrated, n = n_reference, U = max(U)),
    B_calibrated$cell_anchor$I_empirical,
    tolerance = 1e-10
  )),
  all(diff(predict(B_calibrated, n = cell_n, U = max(U))) >= -1e-10)
)
