# Complete q calibration and deterministic evaluation for the Gaussian spike model.

scscale_project_reference_spikes <- function(matrix, reference_vectors, tau2) {
  reference_vectors <- as_dense_matrix(reference_vectors)
  if (nrow(matrix) != nrow(reference_vectors)) {
    stop(
      "matrix and reference_vectors must contain the same ordered features.",
      call. = FALSE
    )
  }
  if (!is.finite(tau2) || tau2 <= 0) {
    stop("tau2 must be positive and finite.", call. = FALSE)
  }

  scores <- as.matrix(crossprod(reference_vectors, matrix))
  center <- attr(matrix, "row_center")
  if (!is.null(center)) {
    projected_center <- as.numeric(crossprod(reference_vectors, center))
    scores <- sweep(scores, 1L, projected_center, "-")
  }
  projected_variance <- rowSums(scores^2) / ncol(matrix)
  excess_variance <- pmax(projected_variance - tau2, 0)
  data.frame(
    component = seq_len(ncol(reference_vectors)),
    projected_variance = projected_variance,
    excess_variance = excess_variance,
    q_X = excess_variance / tau2
  )
}

scscale_fixed_residual_tau2 <- function(matrix, reference_vectors) {
  reference_vectors <- as_dense_matrix(reference_vectors)
  p <- nrow(matrix)
  n <- ncol(matrix)
  r <- ncol(reference_vectors)
  if (nrow(reference_vectors) != p) {
    stop(
      "matrix and reference_vectors must contain the same ordered features.",
      call. = FALSE
    )
  }
  if (r < 1L || r >= p) {
    stop("reference_vectors must contain between 1 and p - 1 directions.", call. = FALSE)
  }

  scores <- as.matrix(crossprod(reference_vectors, matrix))
  center <- attr(matrix, "row_center")
  centered_ss <- sum(matrix^2)
  if (!is.null(center)) {
    center <- as.numeric(center)
    scores <- sweep(
      scores,
      1L,
      as.numeric(crossprod(reference_vectors, center)),
      "-"
    )
    centered_ss <- centered_ss - n * sum(center^2)
  }
  residual_ss <- max(centered_ss - sum(scores^2), 0)
  tau2 <- residual_ss / (n * (p - r))
  if (!is.finite(tau2) || tau2 <= 0) {
    stop("The fixed residual noise estimate is not positive and finite.", call. = FALSE)
  }
  tau2
}

scscale_fit_umi_q_parameters <- function(calibration, reference_U) {
  components <- sort(unique(calibration$component))
  rows <- lapply(components, function(component) {
    df <- calibration[calibration$component == component, , drop = FALSE]
    keep <- is.finite(df$U) & is.finite(df$q_X)
    df <- df[keep, , drop = FALSE]
    if (nrow(df) < 2L) {
      stop(
        "At least two finite UMI calibration points are required per component.",
        call. = FALSE
      )
    }
    model <- stats::lm(q_X ~ 0 + U, data = df)
    coefficients <- stats::coef(model)
    intercept <- 0
    slope <- unname(coefficients[["U"]])
    if (!is.finite(slope)) slope <- 0
    slope <- pmax(slope, 0)
    r_squared <- summary(model)$r.squared
    q_fit_reference <- slope * reference_U
    data.frame(
      component = component,
      q_intercept_a = intercept,
      q_slope_b = slope,
      q_r_squared = r_squared,
      q_fit_reference = q_fit_reference,
      n_depth_observations = nrow(df),
      umi_intercept_a = intercept,
      umi_slope_b = slope,
      umi_r_squared = r_squared,
      n_umi_observations = nrow(df)
    )
  })
  do.call(rbind, rows)
}

scscale_spike_fit <- function(
  x,
  y,
  umi_sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  umi_replicates = 1,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  mp_max_iter = 300,
  mp_grid_n = 3000,
  x_rank = NULL,
  y_rank = NULL,
  x_force_spikes = FALSE,
  y_force_spikes = FALSE,
  alignment_theta_floor = 1e-6,
  alignment_max_singular = 1 - 1e-8,
  seed = 1,
  use_irlba = TRUE,
  ...
) {
  matched_call <- match.call()
  count_transform <- match.arg(count_transform)
  umi_sampling_rates <- sort(unique(as.numeric(umi_sampling_rates)))
  if (length(umi_sampling_rates) < 2L ||
      any(!is.finite(umi_sampling_rates) |
          umi_sampling_rates <= 0 |
          umi_sampling_rates > 1)) {
    stop(
      "umi_sampling_rates must contain at least two unique values in (0, 1].",
      call. = FALSE
    )
  }
  if (!any(abs(umi_sampling_rates - 1) < sqrt(.Machine$double.eps))) {
    umi_sampling_rates <- sort(c(umi_sampling_rates, 1))
  }
  umi_replicates <- as.integer(umi_replicates)
  if (!is.finite(umi_replicates) || umi_replicates < 1L) {
    stop("umi_replicates must be a positive integer.", call. = FALSE)
  }

  x <- counts_matrix(x)
  y <- counts_matrix(y)
  if (is.null(rownames(x))) rownames(x) <- paste0("feature_", seq_len(nrow(x)))
  if (is.null(rownames(y))) rownames(y) <- paste0("target_", seq_len(nrow(y)))
  if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
  if (is.null(colnames(y))) {
    if (ncol(y) != ncol(x)) {
      stop("y must share cell names with x or have the same number of columns.", call. = FALSE)
    }
    colnames(y) <- colnames(x)
  }
  common_cells <- intersect(colnames(x), colnames(y))
  if (length(common_cells) < 2L) {
    stop("x and y must share at least two cells.", call. = FALSE)
  }
  x <- x[, common_cells, drop = FALSE]
  y <- y[, common_cells, drop = FALSE]

  fit_args <- list(
    input = "counts",
    count_transform = count_transform,
    fit_umi = FALSE,
    mp_max_iter = mp_max_iter,
    mp_grid_n = mp_grid_n,
    store_matrix = TRUE,
    use_irlba = use_irlba
  )
  extra_args <- list(...)
  x_fit <- do.call(
    scscale_fit,
    c(list(x = x, r = x_rank, force_spikes = x_force_spikes), fit_args, extra_args)
  )
  y_fit <- do.call(
    scscale_fit,
    c(list(x = y, r = y_rank, force_spikes = y_force_spikes), fit_args, extra_args)
  )
  r_X <- sum(x_fit$spikes$is_spike)
  r_Y <- sum(y_fit$spikes$is_spike)
  if (r_X < 1L || r_Y < 1L) {
    stop(
      "The reference fit must detect at least one spike in both x and y.",
      call. = FALSE
    )
  }

  z_X <- right_singular_vectors(x_fit$matrix, r = r_X, use_irlba = use_irlba)
  z_Y <- right_singular_vectors(y_fit$matrix, r = r_Y, use_irlba = use_irlba)
  P <- scscale_subspace_overlap_matrix(z_X, z_Y)
  rownames(P) <- paste0("X", seq_len(nrow(P)))
  colnames(P) <- paste0("Y", seq_len(ncol(P)))
  reference_vectors <- x_fit$left_vectors[, seq_len(r_X), drop = FALSE]

  x_calibration_counts <- x[x_fit$features, , drop = FALSE]
  reference_U <- stats::median(colSums(x_calibration_counts))
  calibration_rows <- vector(
    "list",
    length(umi_sampling_rates) * umi_replicates
  )
  rate_fit_summaries <- vector("list", length(calibration_rows))
  row_i <- 0L
  for (rate_i in seq_along(umi_sampling_rates)) {
    rho <- umi_sampling_rates[[rate_i]]
    for (replicate in seq_len(umi_replicates)) {
      row_i <- row_i + 1L
      rate_seed <- seed + (rate_i - 1L) * umi_replicates + replicate
      x_rate <- if (abs(rho - 1) < sqrt(.Machine$double.eps)) {
        x_calibration_counts
      } else {
        scscale_downsample_counts_fraction(
          x_calibration_counts,
          fraction = rho,
          seed = rate_seed
        )
      }
      rate_matrix <- if (abs(rho - 1) < sqrt(.Machine$double.eps)) {
        x_fit$matrix
      } else {
        scscale_prepare_matrix(
          x_rate,
          input = "counts",
          target_depth = x_fit$normalization$target_depth %||% 1e4,
          count_transform = count_transform,
          center = x_fit$normalization$center %||% TRUE,
          scale = x_fit$normalization$scale %||% FALSE,
          n_features = NULL
        )$matrix
      }
      tau2_fixed <- scscale_fixed_residual_tau2(
        rate_matrix,
        reference_vectors = reference_vectors
      )
      projected <- scscale_project_reference_spikes(
        rate_matrix,
        reference_vectors = reference_vectors,
        tau2 = tau2_fixed
      )
      U <- stats::median(colSums(x_rate))
      calibration_rows[[row_i]] <- cbind(
        data.frame(
          sampling_rate = rho,
          U = U,
          replicate = replicate,
          tau2 = tau2_fixed
        ),
        projected
      )
      rate_fit_summaries[[row_i]] <- data.frame(
        sampling_rate = rho,
        U = U,
        replicate = replicate,
        total_umi_observed = sum(x_rate),
        median_umi_observed = stats::median(colSums(x_rate)),
        reference_rank = r_X,
        tau2 = tau2_fixed
      )
    }
  }
  calibration <- do.call(rbind, calibration_rows)
  rate_fit_summaries <- do.call(rbind, rate_fit_summaries)
  umi_parameters <- scscale_fit_umi_q_parameters(
    calibration,
    reference_U = reference_U
  )

  x_parameters <- data.frame(
    component = seq_len(r_X),
    lambda = x_fit$spikes$lambda[seq_len(r_X)],
    tau2 = x_fit$spikes$tau2[seq_len(r_X)],
    q_ref = x_fit$spikes$q_X[seq_len(r_X)],
    theta_ref = x_fit$theta_X[seq_len(r_X)],
    theta_infinity_ref = x_fit$theta_infinity[seq_len(r_X)]
  )
  x_parameters <- merge(
    x_parameters,
    umi_parameters,
    by = "component",
    all.x = TRUE,
    sort = FALSE
  )
  x_parameters <- x_parameters[order(x_parameters$component), , drop = FALSE]
  y_parameters <- data.frame(
    component = seq_len(r_Y),
    lambda = y_fit$spikes$lambda[seq_len(r_Y)],
    tau2 = y_fit$spikes$tau2[seq_len(r_Y)],
    q_ref = y_fit$spikes$q_X[seq_len(r_Y)],
    theta_ref = y_fit$theta_X[seq_len(r_Y)],
    theta_infinity_ref = y_fit$theta_infinity[seq_len(r_Y)]
  )

  theta_X_alignment_ref <- scscale_recoverability(
    x_parameters$q_fit_reference,
    c_X = x_fit$c_X
  )
  theta_Y_alignment_ref <- y_parameters$theta_ref
  calibrated_alignment <- scscale_calibrate_reference_alignment(
    P,
    theta_X = theta_X_alignment_ref,
    theta_Y = theta_Y_alignment_ref,
    theta_floor = alignment_theta_floor,
    max_singular = alignment_max_singular
  )
  P_sigma <- svd(P, nu = 0, nv = 0)$d
  P_gamma <- pmin(pmax(P_sigma^2, 0), 1 - 1e-12)
  reference_empirical_mi <- -0.5 * sum(log1p(-P_gamma))

  out <- list(
    call = matched_call,
    reference = list(
      n = length(common_cells),
      p_X = x_fit$p,
      p_Y = y_fit$p,
      c_X = x_fit$c_X,
      c_Y = y_fit$c_X,
      total_umi_X = sum(x_calibration_counts),
      median_umi_per_cell_X = reference_U
    ),
    x = list(
      rank = r_X,
      parameters = x_parameters,
      fit = x_fit,
      z = z_X
    ),
    y = list(
      rank = r_Y,
      parameters = y_parameters,
      fit = y_fit,
      z = z_Y
    ),
    alignment = list(
      P = P,
      K = calibrated_alignment$K,
      K_raw = calibrated_alignment$K_raw,
      r_X = r_X,
      r_Y = r_Y,
      theta_X_reference = theta_X_alignment_ref,
      theta_Y_reference = theta_Y_alignment_ref,
      theta_floor = calibrated_alignment$theta_floor,
      max_singular = calibrated_alignment$max_singular,
      singular_raw = calibrated_alignment$singular_raw,
      singular_regularized = calibrated_alignment$singular_regularized,
      singular_clipped = calibrated_alignment$singular_clipped,
      reference_reconstruction_error = calibrated_alignment$reference_reconstruction_error,
      reference_empirical_mi = reference_empirical_mi,
      definition = "K = Theta_X,ref^(-1/2) P_obs Theta_Y,ref^(-1/2), with singular-value clipping"
    ),
    umi = list(
      sampling_rates = umi_sampling_rates,
      U = umi_sampling_rates * reference_U,
      reference_U = reference_U,
      replicates = umi_replicates,
      parameters = umi_parameters,
      formula = "q_X,k(U) = b_k * U",
      fit_method = "ordinary least squares through the origin",
      anchored_at_zero = TRUE,
      anchored_at_full_depth = FALSE,
      component_tracking = "fixed full-depth reference feature directions",
      noise_estimator = "fixed full-depth residual subspace"
    ),
    diagnostics = list(
      umi_calibration = calibration,
      umi_rate_fits = rate_fit_summaries,
      mp_history_X = x_fit$bulk$history,
      mp_history_Y = y_fit$bulk$history
    ),
    preprocessing = list(
      count_transform = count_transform,
      features_X = x_fit$features,
      features_Y = y_fit$features,
      cells = common_cells
    )
  )
  class(out) <- "scscale_spike_fit"
  ceiling <- scscale_spike_infinity(out)
  out$I_infinity <- ceiling$mi
  out$ceiling <- list(
    theta_X = ceiling$theta_X,
    theta_Y = ceiling$theta_Y,
    gamma = ceiling$gamma,
    I_infinity = ceiling$mi,
    definition = "theta_X,k = 1 with fitted reference theta_Y,k and calibrated latent alignment K"
  )
  out$reference_state <- scscale_spike_state(
    out,
    n = out$reference$n,
    sampling_rate = 1
  )
  out$alignment$reference_model_mi <- out$reference_state$I_theory
  out$alignment$reference_mi_difference <- (
    out$reference_state$I_theory - out$alignment$reference_empirical_mi
  )
  out
}

scscale_calibrate_reference_alignment <- function(
  P,
  theta_X,
  theta_Y,
  theta_floor = 1e-6,
  max_singular = 1 - 1e-8
) {
  P <- as_dense_matrix(P)
  r_X <- min(nrow(P), length(theta_X))
  r_Y <- min(ncol(P), length(theta_Y))
  if (r_X < 1L || r_Y < 1L) {
    stop("P and both theta vectors must be non-empty.", call. = FALSE)
  }
  if (!is.finite(theta_floor) || theta_floor <= 0 || theta_floor > 1) {
    stop("theta_floor must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(max_singular) || max_singular <= 0 || max_singular > 1) {
    stop("max_singular must be in (0, 1].", call. = FALSE)
  }

  P_use <- P[seq_len(r_X), seq_len(r_Y), drop = FALSE]
  theta_X_use <- pmin(pmax(theta_X[seq_len(r_X)], theta_floor), 1)
  theta_Y_use <- pmin(pmax(theta_Y[seq_len(r_Y)], theta_floor), 1)
  K_raw <- sweep(P_use, 1L, sqrt(theta_X_use), "/")
  K_raw <- sweep(K_raw, 2L, sqrt(theta_Y_use), "/")

  decomposition <- svd(K_raw)
  singular_regularized <- pmin(decomposition$d, max_singular)
  K <- sweep(
    decomposition$u,
    2L,
    singular_regularized,
    "*"
  ) %*% t(decomposition$v)
  dimnames(K) <- dimnames(P_use)
  dimnames(K_raw) <- dimnames(P_use)

  P_reconstructed <- sweep(K, 1L, sqrt(theta_X_use), "*")
  P_reconstructed <- sweep(P_reconstructed, 2L, sqrt(theta_Y_use), "*")

  list(
    K = K,
    K_raw = K_raw,
    P_reconstructed = P_reconstructed,
    theta_X = theta_X_use,
    theta_Y = theta_Y_use,
    theta_floor = theta_floor,
    max_singular = max_singular,
    singular_raw = decomposition$d,
    singular_regularized = singular_regularized,
    singular_clipped = decomposition$d > max_singular,
    reference_reconstruction_error = sqrt(sum((P_reconstructed - P_use)^2))
  )
}

scscale_spike_alignment_matrix <- function(fit) {
  fit$alignment$K %||% fit$alignment$P
}

scscale_check_spike_fit <- function(fit) {
  if (!inherits(fit, "scscale_spike_fit")) {
    stop("fit must be a scscale_spike_fit object.", call. = FALSE)
  }
  invisible(fit)
}

scscale_spike_infinity <- function(fit, eps = 1e-12) {
  scscale_check_spike_fit(fit)
  K <- scscale_spike_alignment_matrix(fit)
  theta_X <- rep(1, nrow(K))
  theta_Y <- fit$y$parameters$theta_ref[seq_len(ncol(K))]
  details <- scscale_low_rank_mi(
    theta_X,
    theta_Y,
    P = K,
    eps = eps
  )
  list(
    theta_X = theta_X,
    theta_Y = theta_Y,
    gamma = details$gamma,
    mi = details$mi
  )
}

scscale_spike_state <- function(fit, n, sampling_rate, eps = 1e-12) {
  scscale_check_spike_fit(fit)
  n <- as.numeric(n)
  sampling_rate <- as.numeric(sampling_rate)
  if (length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("n must be one positive finite cell number.", call. = FALSE)
  }
  if (length(sampling_rate) != 1L ||
      !is.finite(sampling_rate) ||
      sampling_rate <= 0) {
    stop("sampling_rate must be one positive finite value.", call. = FALSE)
  }

  x_parameters <- fit$x$parameters
  y_parameters <- fit$y$parameters
  U <- sampling_rate * fit$reference$median_umi_per_cell_X
  q_X <- pmax(x_parameters$q_slope_b * U, 0)
  c_X <- fit$reference$p_X / n
  c_Y <- fit$reference$p_Y / n
  q_Y <- y_parameters$q_ref
  theta_X <- scscale_recoverability(q_X, c_X = c_X)
  theta_Y <- y_parameters$theta_ref
  ceiling <- scscale_spike_infinity(fit, eps = eps)
  theta_X_infinity <- ceiling$theta_X
  theta_Y_infinity <- ceiling$theta_Y
  details <- scscale_low_rank_mi(
    theta_X,
    theta_Y,
    P = scscale_spike_alignment_matrix(fit),
    eps = eps
  )

  list(
    n = n,
    sampling_rate = sampling_rate,
    U = U,
    c_X = c_X,
    c_Y = c_Y,
    q_X = q_X,
    q_Y = q_Y,
    theta_X = theta_X,
    theta_Y = theta_Y,
    theta_X_infinity = theta_X_infinity,
    theta_Y_infinity = theta_Y_infinity,
    P = fit$alignment$P,
    K = scscale_spike_alignment_matrix(fit),
    gamma = details$gamma,
    gamma_infinity = ceiling$gamma,
    I_theory = details$mi,
    I_infinity = ceiling$mi,
    effective_rank = length(details$gamma)
  )
}

scscale_spike_grid <- function(fit, n, sampling_rate, eps = 1e-12) {
  scscale_check_spike_fit(fit)
  grid <- expand.grid(
    n = sort(unique(as.numeric(n))),
    sampling_rate = sort(unique(as.numeric(sampling_rate)))
  )
  if (!nrow(grid) ||
      any(!is.finite(grid$n) | grid$n <= 0) ||
      any(!is.finite(grid$sampling_rate) | grid$sampling_rate <= 0)) {
    stop("n and sampling_rate must contain positive finite values.", call. = FALSE)
  }
  states <- lapply(seq_len(nrow(grid)), function(i) {
    scscale_spike_state(
      fit,
      n = grid$n[[i]],
      sampling_rate = grid$sampling_rate[[i]],
      eps = eps
    )
  })
  grid$U <- vapply(states, `[[`, numeric(1), "U")
  grid$c_X <- vapply(states, `[[`, numeric(1), "c_X")
  grid$c_Y <- vapply(states, `[[`, numeric(1), "c_Y")
  grid$I_theory <- vapply(states, `[[`, numeric(1), "I_theory")
  grid$I_infinity <- vapply(states, `[[`, numeric(1), "I_infinity")
  grid$effective_rank <- vapply(states, `[[`, integer(1), "effective_rank")
  attr(grid, "states") <- states
  grid
}

scscale_cell_scaling <- function(fit = NULL, n = NULL, ..., q_X = NULL) {
  if (inherits(fit, "scscale_spike_fit")) {
    if (is.null(n)) stop("n is required.", call. = FALSE)
    out <- scscale_spike_grid(fit, n = n, sampling_rate = 1)
    return(out)
  }

  legacy_input <- fit %||% q_X
  dots <- list(...)
  if (!is.null(n) && is.null(dots$n_grid)) dots$n_grid <- n
  do.call(
    scscale_cell_scaling_legacy,
    c(list(q_X = legacy_input), dots)
  )
}

scscale_umi_scaling <- function(
  fit = NULL,
  sampling_rate = NULL,
  U = NULL,
  ...,
  counts = NULL
) {
  if (inherits(fit, "scscale_spike_fit")) {
    if (!is.null(sampling_rate) && !is.null(U)) {
      stop("Provide only one of sampling_rate or U.", call. = FALSE)
    }
    if (is.null(sampling_rate) && is.null(U)) {
      sampling_rate <- fit$umi$sampling_rates
    } else if (!is.null(U)) {
      sampling_rate <- as.numeric(U) / fit$reference$median_umi_per_cell_X
    }
    out <- scscale_spike_grid(
      fit,
      n = fit$reference$n,
      sampling_rate = sampling_rate
    )
    out$U <- out$sampling_rate * fit$reference$median_umi_per_cell_X
    return(out)
  }

  legacy_counts <- fit %||% counts
  dots <- list(...)
  if (!is.null(sampling_rate) && is.null(dots$sampling_rates)) {
    dots$sampling_rates <- sampling_rate
  }
  if (!is.null(U) && is.null(dots$U_grid)) dots$U_grid <- U
  do.call(
    scscale_umi_scaling_legacy,
    c(list(counts = legacy_counts), dots)
  )
}

scscale_joint_scaling <- function(
  fit,
  n,
  sampling_rate = NULL,
  U = NULL,
  eps = 1e-12
) {
  scscale_check_spike_fit(fit)
  if (!is.null(sampling_rate) && !is.null(U)) {
    stop("Provide only one of sampling_rate or U.", call. = FALSE)
  }
  if (is.null(sampling_rate) && is.null(U)) {
    sampling_rate <- fit$umi$sampling_rates
  } else if (!is.null(U)) {
    sampling_rate <- as.numeric(U) / fit$reference$median_umi_per_cell_X
  }
  out <- scscale_spike_grid(
    fit,
    n = n,
    sampling_rate = sampling_rate,
    eps = eps
  )
  out$U <- out$sampling_rate * fit$reference$median_umi_per_cell_X
  out
}

print.scscale_spike_fit <- function(x, ...) {
  cat("scScale calibrated spike model\n")
  cat("  cells: ", x$reference$n, "\n", sep = "")
  cat("  X features: ", x$reference$p_X, "  rank: ", x$x$rank, "\n", sep = "")
  cat("  Y features: ", x$reference$p_Y, "  rank: ", x$y$rank, "\n", sep = "")
  cat(
    "  UMI rates: ",
    paste(signif(x$umi$sampling_rates, 3), collapse = ", "),
    "\n",
    sep = ""
  )
  cat("  reference MI: ", signif(x$reference_state$I_theory, 5), "\n", sep = "")
  cat("  shared I_infinity: ", signif(x$I_infinity, 5), "\n", sep = "")
  invisible(x)
}
