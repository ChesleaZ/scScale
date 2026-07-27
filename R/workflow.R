# Paired-modality helpers for spectral MI scaling.

scscale_is_one_hot_like <- function(y, tol = 1e-8) {
  y <- as_dense_matrix(y)
  if (!all(is.finite(y))) return(FALSE)
  vals <- unique(as.numeric(y))
  vals <- vals[is.finite(vals)]
  if (!all(abs(vals) < tol | abs(vals - 1) < tol)) return(FALSE)
  col_sums <- colSums(y)
  all(abs(col_sums - 1) < tol)
}

scscale_matrix_rank <- function(y, tol = 1e-8) {
  y <- as_dense_matrix(y)
  if (!length(y)) return(0L)
  qr(y, tol = tol)$rank
}

scscale_target_right_singular_vectors <- function(y, r, use_irlba = TRUE) {
  y <- as_dense_matrix(y)
  r_eff <- min(as.integer(r), nrow(y), ncol(y) - 1L)
  if (!is.finite(r_eff) || r_eff < 1L) return(NULL)
  if (isTRUE(use_irlba) && requireNamespace("irlba", quietly = TRUE)) {
    fit <- irlba::irlba(y, nv = r_eff, nu = 0, maxit = 2000, work = max(3L * r_eff, 30L))
    z <- fit$v[, seq_len(r_eff), drop = FALSE]
  } else {
    fit <- svd(y, nu = 0, nv = r_eff)
    z <- fit$v[, seq_len(r_eff), drop = FALSE]
  }
  rownames(z) <- colnames(y)
  colnames(z) <- paste0("z", seq_len(ncol(z)))
  z
}

scscale_project_fit_to_target_subspace <- function(fit, z_Y, r_X, eps = 1e-12) {
  if (is.null(fit$matrix) || is.null(fit$left_vectors)) return(NULL)
  x <- fit$matrix
  if (!is.null(rownames(z_Y)) && !is.null(colnames(x))) {
    missing <- setdiff(rownames(z_Y), colnames(x))
    if (length(missing)) stop("target subspace contains cells absent from x.", call. = FALSE)
    z_Y <- z_Y[colnames(x), , drop = FALSE]
  } else if (nrow(z_Y) != ncol(x)) {
    stop("target subspace row count must match x cells.", call. = FALSE)
  }

  r_use <- min(as.integer(r_X), ncol(fit$left_vectors), nrow(fit$spikes), nrow(fit$spectrum))
  if (!is.finite(r_use) || r_use < 1L) return(NULL)
  xz <- as.matrix(x %*% z_Y)
  center <- attr(x, "row_center")
  if (!is.null(center)) {
    xz <- xz - tcrossprod(as.numeric(center), colSums(z_Y))
  }
  lambda <- fit$spectrum$lambda[seq_len(r_use)]
  d <- sqrt(pmax(lambda, 0) * fit$n)
  ok <- is.finite(d) & d > eps
  P <- matrix(0, nrow = r_use, ncol = ncol(z_Y))
  P[ok, ] <- sweep(crossprod(fit$left_vectors[, seq_len(r_use), drop = FALSE][, ok, drop = FALSE], xz), 1L, d[ok], "/")
  rownames(P) <- paste0("zX", seq_len(nrow(P)))
  colnames(P) <- colnames(z_Y)
  P
}

scscale_overlap_mi_from_P <- function(P, eps = 1e-12) {
  sigma <- svd(as_dense_matrix(P), nu = 0, nv = 0)$d
  gamma <- pmin(pmax(sigma^2, 0), 1 - eps)
  list(mi = -0.5 * sum(log1p(-gamma)), gamma = gamma, r_eff = length(gamma))
}

scscale_detect_target_mode <- function(
  y,
  low_rank_max = 50,
  rank_tol = 1e-8
) {
  if (!(is.matrix(y) || inherits(y, "Matrix"))) return("subspace")
  y <- counts_matrix(y)
  p <- nrow(y)
  n <- ncol(y)
  if (p < 1L || n < 2L) return("spike")
  if (scscale_is_one_hot_like(y, tol = rank_tol)) return("subspace")
  rank_y <- scscale_matrix_rank(y, tol = rank_tol)
  full_rank <- min(p, n)
  if (rank_y <= low_rank_max && rank_y < full_rank) return("subspace")
  "spike"
}

scscale_target_fit <- function(
  y,
  input = c("matrix", "labels", "numeric"),
  rank = NULL,
  center = TRUE,
  scale = FALSE,
  drop_empty = TRUE,
  rank_tol = 1e-8,
  use_irlba = TRUE
) {
  input <- match.arg(input)
  if (input == "labels") {
    Y <- label_target_matrix(y, center = center, drop_empty = drop_empty)
  } else if (input == "numeric") {
    Y <- numeric_target_matrix(y, center = center, scale = scale)
  } else {
    Y <- counts_matrix(y)
    if (is.null(colnames(Y))) colnames(Y) <- paste0("cell_", seq_len(ncol(Y)))
    Y <- as_dense_matrix(Y)
    if (isTRUE(center)) Y <- sweep(Y, 1L, rowMeans(Y), "-")
    if (isTRUE(scale)) {
      row_sd <- apply(Y, 1L, stats::sd)
      row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
      Y <- sweep(Y, 1L, row_sd, "/")
    }
    Y[!is.finite(Y)] <- 0
  }
  if (is.null(colnames(Y))) colnames(Y) <- paste0("cell_", seq_len(ncol(Y)))
  r_rank <- scscale_matrix_rank(Y, tol = rank_tol)
  r_use <- rank %||% r_rank
  r_use <- min(as.integer(r_use), r_rank, nrow(Y), ncol(Y) - 1L)
  if (!is.finite(r_use) || r_use < 1L) {
    stop("target subspace has no positive rank after preprocessing.", call. = FALSE)
  }
  z <- scscale_target_right_singular_vectors(Y, r = r_use, use_irlba = use_irlba)
  if (is.null(z) || ncol(z) < 1L) {
    stop("could not compute target right singular vectors.", call. = FALSE)
  }
  r_use <- ncol(z)
  object <- list(
    input = input,
    mode = "subspace",
    n = ncol(Y),
    p = nrow(Y),
    c_X = nrow(Y) / ncol(Y),
    c_ratio = nrow(Y) / ncol(Y),
    matrix = Y,
    z = z,
    r = r_use,
    rank = r_rank,
    theta_X = rep(1, r_use),
    theta = rep(1, r_use),
    theta_infinity = rep(1, r_use),
    spikes = data.frame(
      rank = seq_len(r_use),
      lambda = NA_real_,
      eigenvalue = NA_real_,
      lambda_tilde = NA_real_,
      is_spike = TRUE,
      d2_X = Inf,
      d2 = Inf,
      d = Inf,
      q_X = Inf,
      q = Inf,
      tau2 = NA_real_,
      c_X = nrow(Y) / ncol(Y),
      c_ratio = nrow(Y) / ncol(Y),
      lambda_plus = NA_real_,
      mp_upper = NA_real_
    ),
    bulk = list(
      method = "deterministic_subspace",
      tau2 = NA_real_,
      lambda_plus = NA_real_,
      mp_upper = NA_real_,
      spike_index = seq_len(r_use),
      bulk_index = integer(0)
    )
  )
  class(object) <- c("scscale_subspace_fit", "scscale_fit")
  object
}

scscale_pair_fit <- function(
  x,
  y,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  mp_max_iter = 300,
  mp_grid_n = 3000,
  y_mode = c("spike", "subspace", "auto"),
  y_rank = NULL,
  y_center = TRUE,
  y_scale = FALSE,
  y_input = c("auto", "matrix", "labels", "numeric"),
  auto_low_rank_max = 50,
  rank_tol = 1e-8,
  x_rank = NULL,
  x_force_spikes = FALSE,
  use_irlba = TRUE,
  ...
) {
  count_transform <- match.arg(count_transform)
  y_mode <- match.arg(y_mode)
  y_input <- match.arg(y_input)

  x_fit <- scscale_fit(
    x,
    count_transform = count_transform,
    fit_umi = FALSE,
    mp_max_iter = mp_max_iter,
    mp_grid_n = mp_grid_n,
    r = x_rank,
    force_spikes = x_force_spikes,
    store_matrix = TRUE,
    use_irlba = use_irlba,
    ...
  )
  if (y_mode == "auto") {
    y_mode <- scscale_detect_target_mode(y, low_rank_max = auto_low_rank_max, rank_tol = rank_tol)
  }
  if (y_mode == "subspace" && y_input == "auto") {
    y_input <- if (is.matrix(y) || inherits(y, "Matrix")) {
      "matrix"
    } else if (is.numeric(y) || is.integer(y)) {
      "numeric"
    } else {
      "labels"
    }
  }
  if (y_mode == "subspace") {
    x_cells <- colnames(x_fit$matrix)
    if (!(is.matrix(y) || inherits(y, "Matrix"))) {
      y <- align_vector_to_cells(y, x_cells, "y")
    } else {
      y <- counts_matrix(y)
      if (is.null(colnames(y))) {
        if (ncol(y) != length(x_cells)) {
          stop("subspace target matrix must share cell names with x or have the same number of columns.", call. = FALSE)
        }
        colnames(y) <- x_cells
      } else {
        missing_cells <- setdiff(x_cells, colnames(y))
        if (length(missing_cells)) {
          stop("subspace target matrix is missing ", length(missing_cells), " x cells.", call. = FALSE)
        }
        y <- y[, x_cells, drop = FALSE]
      }
    }
  }
  y_fit <- if (y_mode == "subspace") {
    scscale_target_fit(
      y,
      input = y_input,
      rank = y_rank,
      center = y_center,
      scale = y_scale,
      rank_tol = rank_tol,
      use_irlba = use_irlba
    )
  } else {
    scscale_fit(
      y,
      count_transform = count_transform,
      fit_umi = FALSE,
      mp_max_iter = mp_max_iter,
      mp_grid_n = mp_grid_n,
      store_matrix = TRUE,
      use_irlba = use_irlba,
      ...
    )
  }

  r_X <- sum(x_fit$spikes$is_spike)
  r_Y <- if (inherits(y_fit, "scscale_subspace_fit")) y_fit$r else sum(y_fit$spikes$is_spike)
  if (inherits(y_fit, "scscale_subspace_fit")) {
    z_X <- NULL
    z_Y <- y_fit$z
    P <- scscale_project_fit_to_target_subspace(x_fit, z_Y, r_X = r_X)
    if (is.null(P)) {
      z_X <- right_singular_vectors(x_fit$matrix, r = r_X, use_irlba = use_irlba)
      P <- scscale_subspace_overlap_matrix(z_X, z_Y)
    }
    mi <- scscale_mi(x_fit, y_fit, P = P, empirical = FALSE, use_irlba = use_irlba)
    empirical_overlap <- scscale_overlap_mi_from_P(P)
    mi$I_empirical <- empirical_overlap$mi
    mi$empirical <- list(
      available = TRUE,
      I_empirical = empirical_overlap$mi,
      mi = empirical_overlap$mi,
      gamma = empirical_overlap$gamma,
      r_eff = empirical_overlap$r_eff,
      method = "left_eigen_projection"
    )
  } else {
    z_X <- right_singular_vectors(x_fit$matrix, r = r_X, use_irlba = use_irlba)
    z_Y <- right_singular_vectors(y_fit$matrix, r = r_Y, use_irlba = use_irlba)
    P <- scscale_subspace_overlap_matrix(z_X, z_Y)
    mi <- scscale_mi(x_fit, y_fit, P = P, empirical = TRUE, use_irlba = use_irlba)
  }

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
    count_transform = count_transform,
    y_mode = y_mode
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

  curve <- do.call(rbind, rows)
  reference_ranks <- seq_len(pair$r_X)
  scaling <- do.call(rbind, lapply(seq_along(sampling_rates), function(i) {
    fit_rho <- fit_by_rate[[i]]
    q_all <- fit_rho$spikes$q_X[reference_ranks]
    theta_all <- fit_rho$theta_X[reference_ranks]
    data.frame(
      sampling_rate = sampling_rates[[i]],
      total_umi_observed = curve$total_umi_observed[[i]],
      rank = reference_ranks,
      q_X = q_all,
      theta_X = theta_all,
      is_detected = reference_ranks <= curve$r_X[[i]],
      c_X = fit_rho$c_X
    )
  }))
  q_fit_rate <- do.call(rbind, lapply(split(scaling, scaling$rank), function(df) {
    can_fit <- nrow(df) >= 2L && length(unique(df$sampling_rate)) >= 2L
    if (can_fit) {
      model <- stats::lm(q_X ~ sampling_rate, data = df)
      intercept <- unname(stats::coef(model)[[1]])
      slope <- unname(stats::coef(model)[[2]])
      r_squared <- summary(model)$r.squared
      if (!is.finite(r_squared)) r_squared <- NA_real_
    } else {
      intercept <- df$q_X[[1]]
      slope <- NA_real_
      r_squared <- NA_real_
    }
    data.frame(
      rank = df$rank[[1]],
      q_rate_intercept = intercept,
      q_rate_slope = slope,
      q_rate_r2 = r_squared,
      n_rates = nrow(df),
      n_detected_rates = sum(df$is_detected)
    )
  }))
  q_fit_rate <- q_fit_rate[order(q_fit_rate$rank), , drop = FALSE]
  scaling <- merge(
    scaling,
    q_fit_rate,
    by = "rank",
    all.x = TRUE,
    sort = FALSE
  )
  scaling$q_X_rate_hat <- with(
    scaling,
    q_rate_intercept + q_rate_slope * sampling_rate
  )
  scaling <- scaling[order(scaling$sampling_rate, scaling$rank), , drop = FALSE]

  q_total <- stats::aggregate(
    scaling["q_X"],
    by = scaling["sampling_rate"],
    FUN = sum
  )
  names(q_total)[names(q_total) == "q_X"] <- "q_total"
  can_fit_total <- nrow(q_total) >= 2L &&
    length(unique(q_total$sampling_rate)) >= 2L
  if (can_fit_total) {
    total_model <- stats::lm(q_total ~ sampling_rate, data = q_total)
    total_intercept <- unname(stats::coef(total_model)[[1]])
    total_slope <- unname(stats::coef(total_model)[[2]])
    total_r_squared <- summary(total_model)$r.squared
    if (!is.finite(total_r_squared)) total_r_squared <- NA_real_
    q_total$q_total_rate_hat <- unname(
      stats::predict(total_model, newdata = q_total)
    )
  } else {
    total_intercept <- q_total$q_total[[1]]
    total_slope <- NA_real_
    total_r_squared <- NA_real_
    q_total$q_total_rate_hat <- NA_real_
  }
  q_total_rate_fit <- data.frame(
    q_total_rate_intercept = total_intercept,
    q_total_rate_slope = total_slope,
    q_total_rate_r2 = total_r_squared,
    n_rates = nrow(q_total)
  )

  out <- list(
    curve = curve,
    q_by_rate = q_by_rate,
    fit_by_rate = fit_by_rate,
    P_by_rate = P_by_rate,
    sampling_rates = sampling_rates,
    scaling = scaling,
    q_fit_rate = q_fit_rate,
    q_total = q_total,
    q_total_rate_fit = q_total_rate_fit
  )
  class(out) <- "scscale_umi_mi"
  out
}

scscale_linear_umi_mi <- function(
  pair,
  U_grid = NULL,
  sampling_rates = NULL,
  reference_U = NULL,
  intercept,
  slope,
  r = pair$r_X,
  P = pair$P,
  clamp_rate = TRUE
) {
  if (!inherits(pair, "scscale_pair_fit")) {
    stop("pair must be a scscale_pair_fit object.", call. = FALSE)
  }
  if (is.null(U_grid) && is.null(sampling_rates)) {
    stop("Provide U_grid or sampling_rates.", call. = FALSE)
  }
  if (is.null(sampling_rates)) {
    factor <- scscale_linear_umi_factor(
      U = U_grid,
      reference_U = reference_U,
      intercept = intercept,
      slope = slope,
      clamp_rate = clamp_rate
    )
    sampling_rates <- if (is.null(reference_U)) rep(NA_real_, length(U_grid)) else U_grid / reference_U
  } else {
    factor <- scscale_linear_umi_factor(
      sampling_rate = sampling_rates,
      intercept = intercept,
      slope = slope,
      clamp_rate = clamp_rate
    )
    if (is.null(U_grid)) {
      U_grid <- if (is.null(reference_U)) sampling_rates else sampling_rates * reference_U
    }
  }

  r_X <- min(as.integer(r), pair$r_X, nrow(pair$x_fit$spikes))
  if (!is.finite(r_X) || r_X < 1L) stop("No pair spike strengths available to scale.", call. = FALSE)
  theta_Y <- pair$y_fit$theta_X[seq_len(pair$r_Y)]
  q_full <- pair$x_fit$spikes$q_X[seq_len(r_X)]
  P_use <- P[seq_len(r_X), seq_len(pair$r_Y), drop = FALSE]

  rows <- vector("list", length(factor))
  q_by_rate <- vector("list", length(factor))
  fit_by_rate <- vector("list", length(factor))
  P_by_rate <- vector("list", length(factor))
  names(q_by_rate) <- names(fit_by_rate) <- names(P_by_rate) <- as.character(sampling_rates)

  for (i in seq_along(factor)) {
    q_X <- q_full * factor[[i]]
    theta_X <- scscale_recoverability(q_X, c_X = pair$x_fit$c_X)
    theta_X_infinity <- scscale_theta_infinity(q_X)
    mi <- scscale_low_rank_mi(theta_X, theta_Y, P = P_use)
    I_infinity <- scscale_low_rank_mi(theta_X_infinity, pair$y_fit$theta_infinity, P = P_use)$mi
    rows[[i]] <- data.frame(
      U = U_grid[[i]],
      sampling_rate = sampling_rates[[i]],
      linear_umi_factor = factor[[i]],
      total_umi_observed = NA_real_,
      r_X = r_X,
      I_theory = mi$I_theory,
      I_infinity = I_infinity,
      tau2 = pair$x_fit$bulk$tau2,
      lambda_plus = pair$x_fit$bulk$lambda_plus,
      selected_iteration = pair$x_fit$bulk$selected_iteration
    )
    q_by_rate[[i]] <- q_X
    fit_by_rate[[i]] <- pair$x_fit
    P_by_rate[[i]] <- P_use
  }

  out <- list(
    curve = do.call(rbind, rows),
    q_by_rate = q_by_rate,
    fit_by_rate = fit_by_rate,
    P_by_rate = P_by_rate,
    sampling_rates = sampling_rates,
    U_grid = U_grid,
    linear_umi_parameter = list(
      intercept = intercept,
      slope = slope,
      reference_U = reference_U,
      clamp_rate = isTRUE(clamp_rate),
      formula = "q_X(U) = q_X(full) * (intercept + slope * rho) / (intercept + slope)"
    )
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
