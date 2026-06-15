#!/usr/bin/env Rscript

# Marchenko-Pastur bulk fitting and BBP spike utilities.

mp_edges <- function(gamma, sigma2 = 1) {
  if (!is.finite(gamma) || !is.finite(sigma2) || gamma <= 0 || sigma2 <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(
    lower = sigma2 * (1 - sqrt(gamma))^2,
    upper = sigma2 * (1 + sqrt(gamma))^2
  )
}

mp_positive_mass <- function(gamma) {
  ifelse(gamma > 1, 1 / gamma, 1)
}

mp_density <- function(x, gamma, sigma2 = 1, positive_only = FALSE) {
  edge <- mp_edges(gamma, sigma2)
  out <- numeric(length(x))
  inside <- x >= edge[["lower"]] & x <= edge[["upper"]] & x > 0
  out[inside] <- sqrt((edge[["upper"]] - x[inside]) * (x[inside] - edge[["lower"]])) /
    (2 * pi * gamma * sigma2 * x[inside])
  if (positive_only) out <- out / mp_positive_mass(gamma)
  out
}

mp_quantiles <- function(
  gamma,
  sigma2,
  probs,
  positive_only = TRUE,
  grid_n = 20000
) {
  edge <- mp_edges(gamma, sigma2)
  if (!all(is.finite(edge))) return(rep(NA_real_, length(probs)))
  x <- seq(max(edge[["lower"]], 0) + 1e-12, edge[["upper"]] - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2, positive_only = positive_only)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  stats::approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

mp_median <- function(gamma, positive_only = TRUE, grid_n = 20000) {
  mp_quantiles(gamma, sigma2 = 1, probs = 0.5, positive_only = positive_only, grid_n = grid_n)
}

fit_mp_median <- function(
  eigenvalues,
  n_cells,
  n_features = length(eigenvalues),
  gamma = n_features / n_cells,
  edge_buffer = 1,
  min_positive = 1e-10,
  positive_only = TRUE,
  grid_n = 20000
) {
  vals <- sort(eigenvalues[is.finite(eigenvalues)], decreasing = TRUE)
  positive <- vals[vals > min_positive]
  if (length(positive) < 10) stop("Need at least 10 positive eigenvalues to fit MP.", call. = FALSE)
  med_theory <- mp_median(gamma, positive_only = positive_only, grid_n = grid_n)
  sigma2 <- stats::median(positive) / med_theory
  edge <- mp_edges(gamma, sigma2)
  threshold <- edge_buffer * edge[["upper"]]
  outliers <- positive[positive > threshold]
  data.frame(
    fit_type = "median_fixed_gamma",
    n_cells = n_cells,
    n_features = n_features,
    n_eigenvalues = length(vals),
    n_positive = length(positive),
    n_fit = length(positive),
    gamma_theory = n_features / n_cells,
    gamma = gamma,
    sigma2 = sigma2,
    mp_lower = edge[["lower"]],
    mp_upper = edge[["upper"]],
    edge_buffer = edge_buffer,
    mp_threshold = threshold,
    n_outliers_above_mp = length(outliers),
    outlier_excess_sum = sum(outliers - edge[["upper"]]),
    lambda1_over_mp_upper = ifelse(length(positive) && edge[["upper"]] > 0, positive[1] / edge[["upper"]], NA_real_),
    objective = NA_real_,
    convergence = NA_integer_,
    iterations = 1L,
    mp_median_theory = med_theory,
    empirical_median = stats::median(positive)
  )
}

fit_mp_quantile <- function(
  eigenvalues,
  gamma = NULL,
  probs = seq(0.05, 0.80, by = 0.025),
  grid_n = 20000,
  min_positive = 1e-10,
  gamma_bounds = c(1e-4, 50)
) {
  vals <- eigenvalues[is.finite(eigenvalues) & eigenvalues > min_positive]
  if (length(vals) < 10) stop("Need at least 10 positive eigenvalues to fit MP.", call. = FALSE)
  empirical_q <- as.numeric(stats::quantile(vals, probs, names = FALSE))

  if (!is.null(gamma)) {
    objective <- function(log_sigma2) {
      sigma2 <- exp(log_sigma2)
      fitted_q <- mp_quantiles(gamma, sigma2, probs, positive_only = TRUE, grid_n = grid_n)
      mean((log(empirical_q) - log(fitted_q))^2)
    }
    starts <- log(unique(pmax(c(stats::median(vals), mean(vals), stats::quantile(vals, 0.75, names = FALSE)), .Machine$double.eps)))
    fits <- lapply(starts, function(start) {
      stats::optim(
        start,
        objective,
        method = "Brent",
        lower = log(min(vals) / 100),
        upper = log(max(vals) * 100)
      )
    })
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
    sigma2 <- unname(exp(best$par))
    edge <- mp_edges(gamma, sigma2)
    return(data.frame(
      fit_type = "fixed_gamma",
      gamma = gamma,
      sigma2 = sigma2,
      mp_lower = edge[["lower"]],
      mp_upper = edge[["upper"]],
      objective = best$value,
      convergence = best$convergence,
      n_fit = length(vals)
    ))
  }

  objective <- function(par) {
    gamma_i <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_quantiles(gamma_i, sigma2, probs, positive_only = TRUE, grid_n = grid_n)
    mean((log(empirical_q) - log(fitted_q))^2)
  }
  starts <- expand.grid(
    gamma = unique(pmax(c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8), gamma_bounds[1])),
    sigma2 = unique(pmax(c(stats::median(vals), mean(vals), stats::quantile(vals, 0.75, names = FALSE)), .Machine$double.eps))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    stats::optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = log(c(gamma_bounds[1], min(vals) / 100)),
      upper = log(c(gamma_bounds[2], max(vals) * 100)),
      control = list(maxit = 1500)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma_hat <- unname(exp(best$par[1]))
  sigma2 <- unname(exp(best$par[2]))
  edge <- mp_edges(gamma_hat, sigma2)
  data.frame(
    fit_type = "free_gamma",
    gamma = gamma_hat,
    sigma2 = sigma2,
    mp_lower = edge[["lower"]],
    mp_upper = edge[["upper"]],
    objective = best$value,
    convergence = best$convergence,
    n_fit = length(vals)
  )
}

fit_mp_bulk_iterative <- function(
  eigenvalues,
  n_cells,
  n_features = length(eigenvalues),
  fit_type = c("free_gamma", "fixed_gamma"),
  edge_buffer = 1,
  max_iter = 5,
  min_positive = 1e-10
) {
  fit_type <- match.arg(fit_type)
  vals <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > min_positive], decreasing = TRUE)
  if (length(vals) < 10) stop("Need at least 10 positive eigenvalues.", call. = FALSE)
  fit_vals <- vals
  fit <- NULL
  gamma_theory <- n_features / n_cells

  for (iter in seq_len(max_iter)) {
    fit <- if (fit_type == "fixed_gamma") {
      fit_mp_quantile(fit_vals, gamma = gamma_theory, min_positive = min_positive)
    } else {
      fit_mp_quantile(fit_vals, gamma = NULL, min_positive = min_positive)
    }
    next_vals <- vals[vals <= edge_buffer * fit$mp_upper]
    if (length(next_vals) < 10 || length(next_vals) == length(fit_vals)) break
    fit_vals <- next_vals
  }

  fit$gamma_theory <- gamma_theory
  fit$n_cells <- n_cells
  fit$n_features <- n_features
  fit$n_positive <- length(vals)
  fit$n_fit <- length(fit_vals)
  fit$n_outliers <- sum(vals > edge_buffer * fit$mp_upper)
  fit$edge_buffer <- edge_buffer
  fit$iterations <- iter
  fit
}

fit_mp_robust <- function(
  eigenvalues,
  n_cells,
  n_features = length(eigenvalues),
  fit_type = c("free_gamma", "fixed_gamma"),
  edge_buffer = 1.01,
  max_iter = 5,
  min_positive = 1e-10
) {
  fit_type <- match.arg(fit_type)
  vals <- sort(eigenvalues[is.finite(eigenvalues)], decreasing = TRUE)
  positive <- vals[vals > min_positive]
  bulk_fit <- fit_mp_bulk_iterative(
    positive,
    n_cells = n_cells,
    n_features = n_features,
    fit_type = fit_type,
    edge_buffer = edge_buffer,
    max_iter = max_iter,
    min_positive = min_positive
  )
  threshold <- edge_buffer * bulk_fit$mp_upper
  outliers <- positive[positive > threshold]
  data.frame(
    fit_type = bulk_fit$fit_type,
    n_cells = n_cells,
    n_features = n_features,
    n_eigenvalues = length(vals),
    n_positive = length(positive),
    n_fit = bulk_fit$n_fit,
    gamma_theory = bulk_fit$gamma_theory,
    gamma = bulk_fit$gamma,
    sigma2 = bulk_fit$sigma2,
    mp_lower = bulk_fit$mp_lower,
    mp_upper = bulk_fit$mp_upper,
    edge_buffer = edge_buffer,
    mp_threshold = threshold,
    n_outliers_above_mp = length(outliers),
    outlier_excess_sum = sum(outliers - bulk_fit$mp_upper),
    lambda1_over_mp_upper = ifelse(length(positive) && is.finite(bulk_fit$mp_upper) && bulk_fit$mp_upper > 0, positive[1] / bulk_fit$mp_upper, NA_real_),
    objective = bulk_fit$objective,
    convergence = bulk_fit$convergence,
    iterations = bulk_fit$iterations
  )
}

bbp_population_to_sample <- function(population_eigenvalue, gamma, sigma2 = 1) {
  theta <- population_eigenvalue / sigma2
  out <- rep(mp_edges(gamma, sigma2)[["upper"]], length(theta))
  detectable <- theta > 1 + sqrt(gamma)
  out[detectable] <- sigma2 * theta[detectable] * (1 + gamma / (theta[detectable] - 1))
  out
}

bbp_sample_to_population <- function(sample_eigenvalue, gamma, sigma2 = 1) {
  y <- sample_eigenvalue / sigma2
  disc <- (y - 1 - gamma)^2 - 4 * gamma
  out <- rep(NA_real_, length(y))
  detectable <- is.finite(disc) & disc >= 0 & sample_eigenvalue > mp_edges(gamma, sigma2)[["upper"]]
  u <- ((y[detectable] - 1 - gamma) + sqrt(disc[detectable])) / 2
  out[detectable] <- sigma2 * (1 + u)
  out
}

fit_bbp_spikes <- function(
  eigenvalues,
  n_cells,
  n_features = length(eigenvalues),
  fit_type = c("median_fixed_gamma", "fixed_gamma", "free_gamma"),
  edge_buffer = 1.01
) {
  fit_type <- match.arg(fit_type)
  vals <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
  bulk_fit <- if (fit_type == "median_fixed_gamma") {
    fit_mp_median(vals, n_cells = n_cells, n_features = n_features, edge_buffer = edge_buffer)
  } else {
    fit_mp_bulk_iterative(
      vals,
      n_cells = n_cells,
      n_features = n_features,
      fit_type = fit_type,
      edge_buffer = edge_buffer
    )
  }
  threshold <- edge_buffer * bulk_fit$mp_upper
  spike_vals <- vals[vals > threshold]
  pop_hat <- bbp_sample_to_population(spike_vals, bulk_fit$gamma, bulk_fit$sigma2)
  spike_strength_d <- sqrt(pmax((pop_hat / bulk_fit$sigma2 - 1) / bulk_fit$gamma, 0))
  spike_df <- data.frame(
    spike_rank = seq_along(spike_vals),
    sample_eigenvalue = spike_vals,
    mp_upper = bulk_fit$mp_upper,
    edge_buffer = edge_buffer,
    sample_over_edge = spike_vals / bulk_fit$mp_upper,
    population_eigenvalue_hat = pop_hat,
    spike_strength_d = spike_strength_d,
    spike_strength_d2 = spike_strength_d^2,
    sigma2 = bulk_fit$sigma2,
    gamma = bulk_fit$gamma
  )
  list(bulk_fit = bulk_fit, spikes = spike_df)
}

mp_curve <- function(gamma, sigma2, n = 1000, positive_only = TRUE) {
  edge <- mp_edges(gamma, sigma2)
  x <- seq(max(edge[["lower"]], 0) + 1e-12, edge[["upper"]] - 1e-12, length.out = n)
  data.frame(
    eigenvalue = x,
    density = mp_density(x, gamma, sigma2, positive_only = positive_only)
  )
}
