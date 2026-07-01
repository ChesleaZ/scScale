# Paired-modality helpers for spectral MI scaling.

scscale_pair_fit <- function(
  x,
  y,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  mp_max_iter = 300,
  mp_grid_n = 3000,
  use_irlba = TRUE,
  ...
) {
  count_transform <- match.arg(count_transform)

  x_fit <- scscale_fit(
    x,
    count_transform = count_transform,
    fit_umi = FALSE,
    mp_max_iter = mp_max_iter,
    mp_grid_n = mp_grid_n,
    store_matrix = TRUE,
    use_irlba = use_irlba,
    ...
  )
  y_fit <- scscale_fit(
    y,
    count_transform = count_transform,
    fit_umi = FALSE,
    mp_max_iter = mp_max_iter,
    mp_grid_n = mp_grid_n,
    store_matrix = TRUE,
    use_irlba = use_irlba,
    ...
  )

  r_X <- sum(x_fit$spikes$is_spike)
  r_Y <- sum(y_fit$spikes$is_spike)
  z_X <- right_singular_vectors(x_fit$matrix, r = r_X, use_irlba = use_irlba)
  z_Y <- right_singular_vectors(y_fit$matrix, r = r_Y, use_irlba = use_irlba)
  P <- scscale_subspace_overlap_matrix(z_X, z_Y)
  mi <- scscale_mi(x_fit, y_fit, r = min(r_X, r_Y), P = P, empirical = TRUE, use_irlba = use_irlba)

  out <- list(
    x_fit = x_fit,
    y_fit = y_fit,
    r_X = r_X,
    r_Y = r_Y,
    z_X = z_X,
    z_Y = z_Y,
    P = P,
    mi = mi,
    I_infinity = scscale_low_rank_mi(x_fit$theta_infinity, y_fit$theta_infinity, P = P)$mi,
    count_transform = count_transform
  )
  class(out) <- "scscale_pair_fit"
  out
}

scscale_umi_mi <- function(
  pair,
  x_counts,
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  count_transform = pair$count_transform %||% "log1p_cpm",
  seed = 1,
  mp_max_iter = 300,
  mp_grid_n = 2000,
  use_irlba = TRUE
) {
  if (!inherits(pair, "scscale_pair_fit")) {
    stop("pair must be a scscale_pair_fit object.", call. = FALSE)
  }
  if (any(!is.finite(sampling_rates) | sampling_rates <= 0 | sampling_rates > 1)) {
    stop("sampling_rates must be in (0, 1].", call. = FALSE)
  }

  rows <- vector("list", length(sampling_rates))
  q_by_rate <- vector("list", length(sampling_rates))
  fit_by_rate <- vector("list", length(sampling_rates))
  P_by_rate <- vector("list", length(sampling_rates))
  names(q_by_rate) <- names(fit_by_rate) <- names(P_by_rate) <- as.character(sampling_rates)

  theta_Y <- pair$y_fit$theta_X[seq_len(pair$r_Y)]
  for (i in seq_along(sampling_rates)) {
    rho <- sampling_rates[i]
    x_sub <- x_counts
    if (rho < 1) {
      x_sub <- scscale_downsample_counts_fraction(x_sub, fraction = rho, seed = seed + i)
    }

    fit_rho <- if (rho == 1) {
      pair$x_fit
    } else {
      scscale_fit(
        x_sub,
        count_transform = count_transform,
        fit_umi = FALSE,
        mp_max_iter = mp_max_iter,
        mp_grid_n = mp_grid_n,
        store_matrix = TRUE,
        use_irlba = use_irlba
      )
    }

    r_X <- sum(fit_rho$spikes$is_spike)
    z_X <- right_singular_vectors(fit_rho$matrix, r = r_X, use_irlba = use_irlba)
    P <- scscale_subspace_overlap_matrix(z_X, pair$z_Y)
    mi <- scscale_low_rank_mi(fit_rho$theta_X[seq_len(r_X)], theta_Y, P = P)
    I_infinity <- scscale_low_rank_mi(fit_rho$theta_infinity, pair$y_fit$theta_infinity, P = P)$mi

    rows[[i]] <- data.frame(
      sampling_rate = rho,
      total_umi_observed = sum(x_sub),
      r_X = r_X,
      I_theory = mi$I_theory,
      I_infinity = I_infinity,
      tau2 = fit_rho$bulk$tau2,
      lambda_plus = fit_rho$bulk$lambda_plus,
      selected_iteration = fit_rho$bulk$selected_iteration
    )
    q_by_rate[[i]] <- fit_rho$spikes$q_X[seq_len(r_X)]
    fit_by_rate[[i]] <- fit_rho
    P_by_rate[[i]] <- P
  }

  out <- list(
    curve = do.call(rbind, rows),
    q_by_rate = q_by_rate,
    fit_by_rate = fit_by_rate,
    P_by_rate = P_by_rate,
    sampling_rates = sampling_rates
  )
  class(out) <- "scscale_umi_mi"
  out
}

scscale_cell_number_mi <- function(
  pair,
  n_grid,
  q_X = pair$x_fit$spikes$q_X[seq_len(pair$r_X)],
  P = pair$P
) {
  if (!inherits(pair, "scscale_pair_fit")) {
    stop("pair must be a scscale_pair_fit object.", call. = FALSE)
  }
  theta_Y <- pair$y_fit$theta_X[seq_len(pair$r_Y)]
  out <- scscale_cell_scaling(q_X, p = pair$x_fit$p, n_grid = n_grid, theta_Y = theta_Y, P = P)
  stats::aggregate(I_theory ~ n + c_X, data = out, FUN = unique)
}

scscale_cell_number_by_umi_mi <- function(
  pair,
  umi,
  n_grid,
  sampling_rates = umi$sampling_rates
) {
  if (!inherits(pair, "scscale_pair_fit")) {
    stop("pair must be a scscale_pair_fit object.", call. = FALSE)
  }
  if (!inherits(umi, "scscale_umi_mi")) {
    stop("umi must be a scscale_umi_mi object.", call. = FALSE)
  }

  rows <- lapply(sampling_rates, function(rho) {
    key <- as.character(rho)
    out <- scscale_cell_number_mi(
      pair,
      n_grid = n_grid,
      q_X = umi$q_by_rate[[key]],
      P = umi$P_by_rate[[key]]
    )
    out$sampling_rate <- rho
    out
  })
  do.call(rbind, rows)
}
