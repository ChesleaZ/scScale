# UMI-anchored effective one-dimensional scaling model.

scscale_effective_signal <- function(q_linear, n, d) {
  q_linear <- as.numeric(q_linear)
  n <- as.numeric(n)
  d <- as.numeric(d)
  if (!length(q_linear) || !length(n) || length(d) != 1L ||
      any(!is.finite(q_linear)) ||
      any(!is.finite(n) | n <= 0) ||
      !is.finite(d) || d <= 0) {
    stop("q_linear must be finite, n and d must be positive, and d must be scalar.", call. = FALSE)
  }
  length_out <- max(length(q_linear), length(n))
  if (length_out %% length(q_linear) != 0L || length_out %% length(n) != 0L) {
    stop("q_linear and n must have compatible lengths.", call. = FALSE)
  }
  q_linear <- rep(q_linear, length.out = length_out)
  n <- rep(n, length.out = length_out)
  q_floor <- sqrt(d / n)
  # Shifted softplus floor. The fixed width q_floor makes the construction
  # parameter-free, monotone, strictly above sqrt(d / n), and asymptotically
  # linear when q_linear is large.
  scaled <- (q_linear - q_floor) / q_floor
  softplus <- pmax(scaled, 0) + log1p(exp(-abs(scaled)))
  q_floor + q_floor * softplus
}

scscale_effective_recoverability <- function(q, n, d) {
  q <- scscale_effective_signal(q, n = n, d = d)
  n <- rep(as.numeric(n), length.out = length(q))
  theta <- (q^2 - d / n) / (q * (q + 1))
  pmin(pmax(theta, 0), 1)
}

scscale_effective_mi <- function(
  n,
  U,
  d,
  I_infinity,
  C,
  theta_Y,
  a,
  b
) {
  parameters <- list(
    d = d,
    I_infinity = I_infinity,
    C = C,
    theta_Y = theta_Y,
    a = a,
    b = b
  )
  if (any(vapply(parameters, length, integer(1)) != 1L)) {
    stop("d and all model parameters must be scalars.", call. = FALSE)
  }
  values <- unlist(parameters, use.names = FALSE)
  if (any(!is.finite(values)) ||
      d <= 0 || C <= 0 || theta_Y <= 0 || theta_Y >= 1 ||
      a < 0 || b <= 0) {
    stop(
      "Require d > 0, C > 0, theta_Y in (0, 1), a >= 0, and b > 0.",
      call. = FALSE
    )
  }
  n <- as.numeric(n)
  U <- as.numeric(U)
  if (!length(n) || any(!is.finite(n) | n <= 0) ||
      !length(U) || any(!is.finite(U) | U < 0)) {
    stop("n must be positive and U must be nonnegative.", call. = FALSE)
  }
  length_out <- max(length(n), length(U))
  if (length_out %% length(n) != 0L || length_out %% length(U) != 0L) {
    stop("n and U must have compatible lengths.", call. = FALSE)
  }
  n <- rep(n, length.out = length_out)
  U <- rep(U, length.out = length_out)
  q <- a + b * U
  theta_X <- scscale_effective_recoverability(q, n = n, d = d)
  I_infinity - C * log((1 - theta_X * theta_Y) / (1 - theta_Y))
}

scscale_effective_profile_linear <- function(
  y, g, w, ceiling_min, ceiling_max = Inf, g_zero = NULL
) {
  C_min <- 1e-12
  if (length(y) != length(g) || length(y) != length(w) ||
      any(!is.finite(y) | !is.finite(g) | !is.finite(w) | w <= 0)) {
    return(NULL)
  }
  if (diff(range(g)) <= sqrt(.Machine$double.eps)) return(NULL)

  w_sum <- sum(w)
  y_bar <- sum(w * y) / w_sum
  g_bar <- sum(w * g) / w_sum
  g_centered_ss <- sum(w * (g - g_bar)^2)
  if (!is.finite(g_centered_ss) || g_centered_ss <= 0) return(NULL)

  C_free <- -sum(w * (g - g_bar) * (y - y_bar)) / g_centered_ss
  I_free <- y_bar + C_free * g_bar

  evaluate <- function(I_infinity, C) {
    fitted <- I_infinity - C * g
    list(
      I_infinity = I_infinity,
      C = C,
      fitted = fitted,
      objective = sum(w * (y - fitted)^2)
    )
  }

  solutions <- list()
  feasible <- function(I_infinity, C) {
    is.finite(I_infinity) && is.finite(C) && C >= C_min &&
      I_infinity >= ceiling_min && I_infinity <= ceiling_max &&
      (is.null(g_zero) || I_infinity - C * g_zero >= -1e-12)
  }
  if (is.finite(C_free) && C_free >= C_min &&
      feasible(I_free, C_free)) {
    solutions[[length(solutions) + 1L]] <- evaluate(I_free, C_free)
  }

  g_ss <- sum(w * g^2)
  C_ceiling <- if (g_ss > 0) {
    sum(w * g * (ceiling_min - y)) / g_ss
  } else {
    0
  }
  C_ceiling <- max(C_ceiling, C_min)
  if (!is.null(g_zero) && g_zero > 0) {
    C_ceiling <- min(C_ceiling, ceiling_min / g_zero)
  }
  if (feasible(ceiling_min, C_ceiling)) {
    solutions[[length(solutions) + 1L]] <- evaluate(ceiling_min, C_ceiling)
  }

  # Active upper boundary I_infinity = ceiling_max.
  if (is.finite(ceiling_max)) {
    C_upper <- if (g_ss > 0) {
      sum(w * g * (ceiling_max - y)) / g_ss
    } else {
      C_min
    }
    C_upper <- max(C_upper, C_min)
    if (!is.null(g_zero) && g_zero > 0) {
      C_upper <- min(C_upper, ceiling_max / g_zero)
    }
    if (feasible(ceiling_max, C_upper)) {
      solutions[[length(solutions) + 1L]] <- evaluate(ceiling_max, C_upper)
    }
  }
  I_at_C_min <- max(ceiling_min, y_bar + C_min * g_bar)
  if (!is.null(g_zero)) I_at_C_min <- max(I_at_C_min, C_min * g_zero)
  if (feasible(I_at_C_min, C_min)) {
    solutions[[length(solutions) + 1L]] <- evaluate(I_at_C_min, C_min)
  }

  # Active physical boundary I(0) = I_infinity - C * g_zero = 0.
  if (!is.null(g_zero) && is.finite(g_zero) && g_zero > 0) {
    h <- g_zero - g
    h_ss <- sum(w * h^2)
    if (h_ss > 0) {
      C_zero <- max(sum(w * h * y) / h_ss, ceiling_min / g_zero, C_min)
      if (is.finite(ceiling_max)) C_zero <- min(C_zero, ceiling_max / g_zero)
      I_zero <- C_zero * g_zero
      if (feasible(I_zero, C_zero)) {
        solutions[[length(solutions) + 1L]] <- evaluate(I_zero, C_zero)
      }
    }
  }

  objectives <- vapply(solutions, function(x) x[["objective"]], numeric(1))
  solutions[[which.min(objectives)]]
}

scscale_effective_summarize_umi <- function(data, U_col, mi_col, weighting) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  missing <- setdiff(c(U_col, mi_col), names(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  U <- as.numeric(data[[U_col]])
  mi <- as.numeric(data[[mi_col]])
  keep <- is.finite(U) & U >= 0 & is.finite(mi)
  U <- U[keep]
  mi <- mi[keep]
  if (!length(U)) stop("No complete UMI and MI observations remain.", call. = FALSE)

  levels <- sort(unique(U))
  rows <- lapply(levels, function(level) {
    values <- mi[U == level]
    data.frame(
      U = level,
      I_mean = mean(values),
      I_variance = if (length(values) > 1L) stats::var(values) / length(values) else NA_real_,
      n_replicates = length(values)
    )
  })
  out <- do.call(rbind, rows)
  if (nrow(out) < 5L) {
    stop("At least five distinct UMI levels are required.", call. = FALSE)
  }

  if (weighting == "inverse_variance" &&
      all(is.finite(out$I_variance)) &&
      any(out$I_variance > 0)) {
    positive <- out$I_variance[out$I_variance > 0]
    variance_floor <- max(stats::median(positive) * 0.05, .Machine$double.eps)
    out$weight <- 1 / pmax(out$I_variance, variance_floor)
    out$weight <- out$weight / mean(out$weight)
    attr(out, "weighting") <- "inverse_variance"
    attr(out, "variance_floor") <- variance_floor
  } else {
    out$weight <- 1
    attr(out, "weighting") <- "equal"
    attr(out, "variance_floor") <- NA_real_
  }
  out
}

scscale_effective_fit <- function(
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
) {
  weighting <- match.arg(weighting)
  d <- as.numeric(d)
  n_reference <- as.numeric(n_reference)
  if (length(d) != 1L || !is.finite(d) || d <= 0 ||
      length(n_reference) != 1L || !is.finite(n_reference) || n_reference <= 0) {
    stop("d and n_reference must be positive finite scalars.", call. = FALSE)
  }
  if (length(n_starts) != 1L || !is.finite(n_starts) ||
      n_starts < 1L || n_starts != as.integer(n_starts)) {
    stop("n_starts must be a positive integer.", call. = FALSE)
  }
  n_starts <- as.integer(n_starts)
  I_infinity_max_factor <- as.numeric(I_infinity_max_factor)
  if (length(I_infinity_max_factor) != 1L ||
      is.na(I_infinity_max_factor) || I_infinity_max_factor < 1) {
    stop("I_infinity_max_factor must be at least 1 or Inf.", call. = FALSE)
  }

  calibration <- scscale_effective_summarize_umi(
    data,
    U_col = U_col,
    mi_col = mi_col,
    weighting = weighting
  )
  weighting_used <- attr(calibration, "weighting")
  variance_floor <- attr(calibration, "variance_floor")
  U <- calibration$U
  y <- calibration$I_mean
  w <- calibration$weight
  U_scale <- max(U)
  if (!is.finite(U_scale) || U_scale <= 0) {
    stop("At least one positive UMI level is required.", call. = FALSE)
  }
  if (is.null(U_reference)) U_reference <- U_scale
  U_reference <- as.numeric(U_reference)
  if (length(U_reference) != 1L ||
      !is.finite(U_reference) || U_reference < 0) {
    stop("U_reference must be one nonnegative finite value.", call. = FALSE)
  }

  c_reference <- d / n_reference
  q_reference_start <- max(1, 1.5 * sqrt(c_reference))
  defaults <- list(
    theta_Y = 0.5,
    a = max(0.1 * q_reference_start, 1e-6),
    b = max(0.9 * q_reference_start / U_scale, 1e-8)
  )
  start <- utils::modifyList(defaults, start %||% list())
  if (!is.finite(start$theta_Y) || start$theta_Y <= 0 || start$theta_Y >= 1 ||
      !is.finite(start$a) || start$a < 0 ||
      !is.finite(start$b) || start$b <= 0) {
    stop("start must have theta_Y in (0, 1), a >= 0, and b > 0.", call. = FALSE)
  }

  pack <- function(theta_Y, a, b) {
    c(
      stats::qlogis(theta_Y),
      log(max(a, 1e-12)),
      log(b * U_scale)
    )
  }
  unpack <- function(par) {
    if (length(par) != 3L || any(!is.finite(par))) return(NULL)
    theta_Y <- stats::plogis(par[[1]])
    a <- exp(pmin(pmax(par[[2]], -30), 30))
    slope_scaled <- exp(pmin(pmax(par[[3]], -30), 30))
    list(theta_Y = theta_Y, a = a, b = slope_scaled / U_scale)
  }

  ceiling_min <- max(y)
  ceiling_max <- I_infinity_max_factor * ceiling_min
  evaluate <- function(par, details = FALSE) {
    shape <- unpack(par)
    if (is.null(shape) || shape$theta_Y <= 1e-8 ||
        shape$theta_Y >= 1 - 1e-8) {
      return(if (details) NULL else .Machine$double.xmax / 100)
    }
    q_linear <- shape$a + shape$b * U
    q <- scscale_effective_signal(q_linear, n = n_reference, d = d)
    theta_X <- scscale_effective_recoverability(q_linear, n = n_reference, d = d)
    g <- log((1 - theta_X * shape$theta_Y) / (1 - shape$theta_Y))
    q_linear_zero <- shape$a
    theta_X_zero <- scscale_effective_recoverability(
      q_linear_zero,
      n = n_reference,
      d = d
    )
    g_zero <- log(
      (1 - theta_X_zero * shape$theta_Y) / (1 - shape$theta_Y)
    )
    linear <- scscale_effective_profile_linear(
      y, g, w, ceiling_min, ceiling_max = ceiling_max, g_zero = g_zero
    )
    if (is.null(linear) || !is.finite(linear$objective)) {
      return(if (details) NULL else .Machine$double.xmax / 100)
    }
    if (!details) return(linear$objective)
    c(shape, list(q_linear = q_linear, q = q, theta_X = theta_X, g = g, linear = linear))
  }

  par0 <- pack(start$theta_Y, start$a, start$b)
  starts <- list(par0)
  if (n_starts > 1L) {
    set.seed(seed)
    for (i in seq_len(n_starts - 1L)) {
      starts[[length(starts) + 1L]] <- par0 + stats::rnorm(3, sd = c(1, 1.5, 1.5))
    }
  }
  optimizations <- lapply(starts, function(initial) {
    stats::optim(initial, evaluate, method = "Nelder-Mead", control = control)
  })
  values <- vapply(optimizations, function(x) x[["value"]], numeric(1))
  best_index <- which.min(values)
  optimization <- optimizations[[best_index]]
  fit_details <- evaluate(optimization$par, details = TRUE)
  if (is.null(fit_details)) {
    stop("Effective scaling optimization did not produce a finite fit.", call. = FALSE)
  }

  coefficients <- c(
    I_infinity = fit_details$linear$I_infinity,
    C = fit_details$linear$C,
    theta_Y = fit_details$theta_Y,
    a = fit_details$a,
    b = fit_details$b
  )
  ceiling_tolerance <- 1e-6 * max(1, ceiling_max)
  I_infinity_boundary <- if (
    is.finite(ceiling_max) &&
      abs(coefficients[["I_infinity"]] - ceiling_max) <= ceiling_tolerance
  ) "upper" else if (
    abs(coefficients[["I_infinity"]] - ceiling_min) <=
      1e-6 * max(1, ceiling_min)
  ) "lower" else "none"
  I_infinity_identified <- I_infinity_boundary != "upper"
  calibration$q_linear <- fit_details$q_linear
  calibration$q <- fit_details$q
  calibration$theta_X <- fit_details$theta_X
  calibration$I_fitted <- fit_details$linear$fitted
  calibration$residual <- calibration$I_mean - calibration$I_fitted
  residual <- calibration$residual
  metrics <- data.frame(
    n_levels = nrow(calibration),
    rss = sum(residual^2),
    weighted_rss = fit_details$linear$objective,
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual))
  )

  out <- list(
    coefficients = coefficients,
    d = d,
    n_reference = n_reference,
    U_reference = U_reference,
    U_scale = U_scale,
    I_infinity_max_factor = I_infinity_max_factor,
    I_infinity_bounds = c(lower = ceiling_min, upper = ceiling_max),
    I_infinity_identified = I_infinity_identified,
    I_infinity_boundary = I_infinity_boundary,
    calibration = calibration,
    weighting = weighting_used %||% weighting,
    variance_floor = variance_floor,
    metrics = metrics,
    optimization = list(
      convergence = optimization$convergence,
      message = optimization$message,
      value = optimization$value,
      counts = optimization$counts,
      n_starts = n_starts,
      selected_start = best_index
    ),
    call = match.call()
  )
  class(out) <- "scscale_effective_fit"
  if (!I_infinity_identified) {
    warning(
      sprintf(
        "I_infinity is not identified over the observed UMI range; the estimate reached the configured ceiling of %g x maximum empirical MI (%g).",
        I_infinity_max_factor,
        ceiling_max
      ),
      call. = FALSE
    )
  }
  out
}

scscale_effective_summarize_cell <- function(
  data,
  n_col,
  U_col,
  mi_col,
  U_reference,
  weighting
) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  missing <- setdiff(c(n_col, mi_col), names(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  n <- as.numeric(data[[n_col]])
  mi <- as.numeric(data[[mi_col]])
  U <- if (U_col %in% names(data)) {
    as.numeric(data[[U_col]])
  } else {
    rep(U_reference, length(n))
  }
  keep <- is.finite(n) & n > 0 & is.finite(U) & U >= 0 & is.finite(mi)
  n <- n[keep]
  U <- U[keep]
  mi <- mi[keep]
  if (!length(n)) stop("No complete cell-number and MI observations remain.", call. = FALSE)

  levels <- sort(unique(n))
  rows <- lapply(levels, function(level) {
    at_level <- n == level
    values <- mi[at_level]
    data.frame(
      n = level,
      U = mean(U[at_level]),
      I_mean = mean(values),
      I_variance = if (length(values) > 1L) stats::var(values) / length(values) else NA_real_,
      n_replicates = length(values)
    )
  })
  out <- do.call(rbind, rows)
  if (nrow(out) < 2L) {
    stop("At least two distinct cell-number levels are required.", call. = FALSE)
  }

  if (weighting == "inverse_variance" &&
      all(is.finite(out$I_variance)) &&
      any(out$I_variance > 0)) {
    positive <- out$I_variance[out$I_variance > 0]
    variance_floor <- max(stats::median(positive) * 0.05, .Machine$double.eps)
    out$weight <- 1 / pmax(out$I_variance, variance_floor)
    out$weight <- out$weight / mean(out$weight)
    attr(out, "weighting") <- "inverse_variance"
    attr(out, "variance_floor") <- variance_floor
  } else {
    out$weight <- 1
    attr(out, "weighting") <- "equal"
    attr(out, "variance_floor") <- NA_real_
  }
  out
}

scscale_effective_A_prediction <- function(object, n, U, gamma) {
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("gamma must be one nonnegative finite value.", call. = FALSE)
  }
  coef <- object$coefficients
  n <- as.numeric(n)
  U <- as.numeric(U)
  length_out <- max(length(n), length(U))
  n <- rep(n, length.out = length_out)
  U <- rep(U, length.out = length_out)
  q_linear <- coef[["a"]] + coef[["b"]] * U
  q <- scscale_effective_signal(q_linear, n = n, d = object$d)
  theta_X <- scscale_effective_recoverability(q_linear, n = n, d = object$d)
  theta_anchor <- scscale_effective_recoverability(
    q_linear,
    n = object$n_reference,
    d = object$d
  )
  mi_from_theta <- function(theta) {
    coef[["I_infinity"]] - coef[["C"]] * log(
      (1 - theta * coef[["theta_Y"]]) / (1 - coef[["theta_Y"]])
    )
  }
  I_theory <- mi_from_theta(theta_X)
  anchor_offset <- if (!is.null(object$cell_anchor)) {
    object$cell_anchor$offset
  } else {
    0
  }
  I_anchor <- mi_from_theta(theta_anchor) + anchor_offset
  I_BBP <- coef[["I_infinity"]] - coef[["C"]] * log(
    1 / (1 - coef[["theta_Y"]])
  )

  n_BBP <- rep(Inf, length_out)
  positive_q <- q > 0
  n_BBP[positive_q] <- object$d / q[positive_q]^2
  rho_infinity <- coef[["theta_Y"]] * q / (q + 1)
  A_n <- rep(NA_real_, length_out)
  valid_A <- positive_q & rho_infinity < 1
  A_n[valid_A] <- coef[["theta_Y"]] * object$d /
    (q[valid_A] * (q[valid_A] + 1) * (1 - rho_infinity[valid_A]))

  I_fitted <- I_theory
  calibrated <- is.finite(A_n) & A_n > 0
  if (any(calibrated)) {
    scaled_A <- gamma * A_n[calibrated]
    I_fitted[calibrated] <- I_anchor[calibrated] + coef[["C"]] * (
      log1p(scaled_A / object$n_reference) -
        log1p(scaled_A / n[calibrated])
    )
  }
  list(
    q_linear = q_linear,
    q = q,
    theta_X = theta_X,
    theta_X_anchor = theta_anchor,
    A_n = A_n,
    n_BBP = n_BBP,
    I_BBP = rep(I_BBP, length_out),
    I_theory = I_theory,
    I_fitted = I_fitted
  )
}

scscale_effective_B_prediction <- function(object, n, U, B_reference) {
  if (length(B_reference) != 1L || !is.finite(B_reference) || B_reference < 0) {
    stop("B_reference must be one nonnegative finite value.", call. = FALSE)
  }
  coef <- object$coefficients
  n <- as.numeric(n)
  U <- as.numeric(U)
  length_out <- max(length(n), length(U))
  n <- rep(n, length.out = length_out)
  U <- rep(U, length.out = length_out)
  q_linear <- coef[["a"]] + coef[["b"]] * U
  q_linear_reference <- coef[["a"]] + coef[["b"]] * object$U_reference
  q <- scscale_effective_signal(
    q_linear,
    n = object$n_reference,
    d = object$d
  )
  q_reference <- scscale_effective_signal(
    q_linear_reference,
    n = object$n_reference,
    d = object$d
  )
  if (!is.finite(q_reference) || q_reference <= 0 || any(!is.finite(q) | q <= 0)) {
    stop("The fitted q(U) must remain positive.", call. = FALSE)
  }
  q_ratio <- q / q_reference
  h <- (1 - coef[["theta_Y"]]) * q_reference
  B_U <- B_reference * (1 + h) /
    (q_ratio * (1 + h * q_ratio))

  I_umi_reference <- scscale_effective_mi(
    n = object$n_reference,
    U = U,
    d = object$d,
    I_infinity = coef[["I_infinity"]],
    C = coef[["C"]],
    theta_Y = coef[["theta_Y"]],
    a = coef[["a"]],
    b = coef[["b"]]
  )
  anchor_offset <- if (!is.null(object$cell_anchor)) {
    object$cell_anchor$offset
  } else {
    0
  }
  I_anchor <- I_umi_reference + anchor_offset
  I_fitted <- I_anchor + coef[["C"]] * (
    log1p(B_U / object$n_reference) - log1p(B_U / n)
  )
  list(
    q_linear = q_linear,
    q = q,
    q_reference = rep(q_reference, length_out),
    q_ratio = q_ratio,
    h = rep(h, length_out),
    B_U = B_U,
    I_anchor = I_anchor,
    I_fitted = I_fitted
  )
}

#' Calibrate the effective cell-number law with a direct reference scale
#'
#' Fits one nonnegative parameter, `B_reference`, at the UMI reference depth.
#' At other UMI depths the scale is transferred by the ratio
#' `A(q(U)) / A(q(U_reference))`, so the arbitrary absolute scale of `A` does
#' not enter the cell-number optimization.
#'
#' @export
scscale_effective_cell_calibrate_B <- function(
  object,
  data,
  n_col = "n",
  U_col = "U",
  mi_col = "I_empirical",
  weighting = c("equal", "inverse_variance"),
  B_bounds = c(1e-12, 1e12)
) {
  if (!inherits(object, "scscale_effective_fit")) {
    stop("object must be a scscale_effective_fit.", call. = FALSE)
  }
  weighting <- match.arg(weighting)
  B_bounds <- as.numeric(B_bounds)
  if (length(B_bounds) != 2L || any(!is.finite(B_bounds)) ||
      B_bounds[[1]] <= 0 || B_bounds[[1]] >= B_bounds[[2]]) {
    stop("B_bounds must contain two increasing positive values.", call. = FALSE)
  }
  calibration <- scscale_effective_summarize_cell(
    data,
    n_col = n_col,
    U_col = U_col,
    mi_col = mi_col,
    U_reference = object$U_reference,
    weighting = weighting
  )
  weighting_used <- attr(calibration, "weighting")
  variance_floor <- attr(calibration, "variance_floor")
  anchor_index <- which(calibration$n == object$n_reference)
  if (length(anchor_index) != 1L) {
    stop(
      "Cell calibration data must contain the UMI reference cell count n_reference.",
      call. = FALSE
    )
  }
  anchor_empirical <- calibration$I_mean[[anchor_index]]
  anchor_model <- scscale_effective_mi(
    n = object$n_reference,
    U = object$U_reference,
    d = object$d,
    I_infinity = object$coefficients[["I_infinity"]],
    C = object$coefficients[["C"]],
    theta_Y = object$coefficients[["theta_Y"]],
    a = object$coefficients[["a"]],
    b = object$coefficients[["b"]]
  )
  object$cell_anchor <- list(
    n = object$n_reference,
    U = object$U_reference,
    I_empirical = anchor_empirical,
    I_model = anchor_model,
    offset = anchor_empirical - anchor_model
  )

  evaluate <- function(B_reference, details = FALSE) {
    prediction <- scscale_effective_B_prediction(
      object,
      n = calibration$n,
      U = calibration$U,
      B_reference = B_reference
    )
    residual <- calibration$I_mean - prediction$I_fitted
    objective <- sum(calibration$weight * residual^2)
    if (!details) return(objective)
    c(prediction, list(residual = residual, objective = objective))
  }
  optimization <- stats::optimize(
    function(log_B) evaluate(exp(log_B)),
    interval = log(B_bounds)
  )
  candidates <- c(0, exp(optimization$minimum))
  objectives <- vapply(candidates, evaluate, numeric(1))
  B_reference <- candidates[[which.min(objectives)]]
  details <- evaluate(B_reference, details = TRUE)

  calibration$q_linear <- details$q_linear
  calibration$q <- details$q
  calibration$q_ratio <- details$q_ratio
  calibration$h <- details$h
  calibration$B_U <- details$B_U
  calibration$I_anchor <- details$I_anchor
  calibration$I_fitted <- details$I_fitted
  calibration$residual <- details$residual
  base_coefficients <- object$coefficients[
    !names(object$coefficients) %in% c("gamma", "B_reference")
  ]
  object$coefficients <- c(base_coefficients, B_reference = B_reference)
  object$cell_mode <- "B_reference"
  object$cell_calibration <- calibration
  object$cell_weighting <- weighting_used %||% weighting
  object$cell_variance_floor <- variance_floor
  object$cell_metrics <- data.frame(
    n_levels = nrow(calibration),
    rss = sum(details$residual^2),
    weighted_rss = details$objective,
    rmse = sqrt(mean(details$residual^2)),
    mae = mean(abs(details$residual))
  )
  object$cell_optimization <- list(
    minimum = B_reference,
    objective = details$objective,
    bounds = B_bounds,
    scale = "B_reference"
  )
  object
}

scscale_effective_cell_calibrate <- function(
  object,
  data,
  n_col = "n",
  U_col = "U",
  mi_col = "I_empirical",
  weighting = c("equal", "inverse_variance"),
  gamma_bounds = c(0, 1e4)
) {
  if (!inherits(object, "scscale_effective_fit")) {
    stop("object must be a scscale_effective_fit.", call. = FALSE)
  }
  weighting <- match.arg(weighting)
  gamma_bounds <- as.numeric(gamma_bounds)
  if (length(gamma_bounds) != 2L || any(!is.finite(gamma_bounds)) ||
      gamma_bounds[[1]] < 0 || gamma_bounds[[1]] >= gamma_bounds[[2]]) {
    stop("gamma_bounds must contain two increasing nonnegative values.", call. = FALSE)
  }

  calibration <- scscale_effective_summarize_cell(
    data,
    n_col = n_col,
    U_col = U_col,
    mi_col = mi_col,
    U_reference = object$U_reference,
    weighting = weighting
  )
  weighting_used <- attr(calibration, "weighting")
  variance_floor <- attr(calibration, "variance_floor")
  anchor_index <- which(calibration$n == object$n_reference)
  if (length(anchor_index) != 1L) {
    stop(
      "Cell calibration data must contain the UMI reference cell count n_reference.",
      call. = FALSE
    )
  }
  anchor_empirical <- calibration$I_mean[[anchor_index]]
  anchor_source <- object
  anchor_source$cell_anchor <- NULL
  anchor_model <- scscale_effective_A_prediction(
    anchor_source,
    n = object$n_reference,
    U = object$U_reference,
    gamma = 0
  )$I_fitted[[1L]]
  object$cell_anchor <- list(
    n = object$n_reference,
    U = object$U_reference,
    I_empirical = anchor_empirical,
    I_model = anchor_model,
    offset = anchor_empirical - anchor_model
  )
  evaluate <- function(gamma, details = FALSE) {
    prediction <- scscale_effective_A_prediction(
      object,
      n = calibration$n,
      U = calibration$U,
      gamma = gamma
    )
    residual <- calibration$I_mean - prediction$I_fitted
    objective <- sum(calibration$weight * residual^2)
    if (!details) return(objective)
    c(prediction, list(residual = residual, objective = objective))
  }

  optimization <- stats::optimize(evaluate, interval = gamma_bounds)
  gamma <- optimization$minimum
  details <- evaluate(gamma, details = TRUE)
  calibration$q <- details$q
  calibration$theta_X <- details$theta_X
  calibration$theta_X_anchor <- details$theta_X_anchor
  calibration$A_n <- details$A_n
  calibration$n_BBP <- details$n_BBP
  calibration$I_BBP <- details$I_BBP
  calibration$I_theory <- details$I_theory
  calibration$I_fitted <- details$I_fitted
  calibration$residual <- details$residual

  base_coefficients <- object$coefficients[names(object$coefficients) != "gamma"]
  object$coefficients <- c(base_coefficients, gamma = gamma)
  object$cell_calibration <- calibration
  object$cell_weighting <- weighting_used %||% weighting
  object$cell_variance_floor <- variance_floor
  object$cell_metrics <- data.frame(
    n_levels = nrow(calibration),
    rss = sum(details$residual^2),
    weighted_rss = details$objective,
    rmse = sqrt(mean(details$residual^2)),
    mae = mean(abs(details$residual))
  )
  object$cell_optimization <- list(
    minimum = gamma,
    objective = optimization$objective,
    bounds = gamma_bounds
  )
  object
}

predict.scscale_effective_fit <- function(
  object,
  newdata = NULL,
  n = NULL,
  U = NULL,
  type = c("response", "all"),
  ...
) {
  type <- match.arg(type)
  if (!inherits(object, "scscale_effective_fit")) {
    stop("object must be a scscale_effective_fit.", call. = FALSE)
  }
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
    missing <- setdiff(c("n", "U"), names(newdata))
    if (length(missing)) {
      stop("newdata must contain n and U columns.", call. = FALSE)
    }
    grid <- data.frame(n = as.numeric(newdata$n), U = as.numeric(newdata$U))
  } else {
    if (is.null(n)) n <- object$n_reference
    if (is.null(U)) U <- object$U_reference
    grid <- expand.grid(
      n = as.numeric(n),
      U = as.numeric(U),
      KEEP.OUT.ATTRS = FALSE
    )
  }
  if (!nrow(grid) ||
      any(!is.finite(grid$n) | grid$n <= 0) ||
      any(!is.finite(grid$U) | grid$U < 0)) {
    stop("Prediction n must be positive and U must be nonnegative.", call. = FALSE)
  }

  coef <- object$coefficients
  if (identical(object$cell_mode, "B_reference")) {
    prediction <- scscale_effective_B_prediction(
      object,
      n = grid$n,
      U = grid$U,
      B_reference = coef[["B_reference"]]
    )
    grid$q_linear <- prediction$q_linear
    grid$q <- prediction$q
    grid$q_ratio <- prediction$q_ratio
    grid$h <- prediction$h
    grid$B_U <- prediction$B_U
    grid$I_anchor <- prediction$I_anchor
    grid$I_fitted <- prediction$I_fitted
    if (type == "response") return(grid$I_fitted)
    return(grid)
  }
  gamma <- if ("gamma" %in% names(coef)) coef[["gamma"]] else 1
  prediction <- scscale_effective_A_prediction(
    object,
    n = grid$n,
    U = grid$U,
    gamma = gamma
  )
  grid$q_linear <- prediction$q_linear
  grid$q <- prediction$q
  grid$theta_X <- prediction$theta_X
  grid$A_n <- prediction$A_n
  grid$n_BBP <- prediction$n_BBP
  grid$I_BBP <- prediction$I_BBP
  grid$I_theory <- prediction$I_theory
  grid$I_fitted <- prediction$I_fitted
  if (type == "response") return(grid$I_fitted)
  grid
}

print.scscale_effective_fit <- function(x, ...) {
  cat("UMI-anchored effective scScale fit\n")
  cat("  representation dimension: ", x$d, "\n", sep = "")
  cat("  reference cells: ", x$n_reference, "\n", sep = "")
  cat("  UMI levels: ", nrow(x$calibration), "\n", sep = "")
  cat("  coefficients:\n")
  print(signif(x$coefficients, 5))
  cat("  RMSE: ", signif(x$metrics$rmse, 5), "\n", sep = "")
  cat("  I_infinity identified: ", x$I_infinity_identified, "\n", sep = "")
  if (!is.null(x$cell_metrics)) {
    cat("  cell-number RMSE: ", signif(x$cell_metrics$rmse, 5), "\n", sep = "")
  }
  invisible(x)
}
