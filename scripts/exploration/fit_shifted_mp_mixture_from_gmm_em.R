library(ggplot2)

gmm_dir <- "outputs/exploration/gmm_spectrum"
init_dir <- "outputs/exploration/shifted_mp_gmm_initialization"
out_dir <- "outputs/exploration/shifted_mp_mixture_from_gmm_em"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_iter <- 200
tol <- 1e-8
density_floor <- 1e-12

mp_density <- function(x, lambda, sigma2) {
  if (!is.finite(lambda) || !is.finite(sigma2) || lambda <= 0 || sigma2 <= 0) {
    return(rep(0, length(x)))
  }
  lower <- sigma2 * (1 - sqrt(lambda))^2
  upper <- sigma2 * (1 + sqrt(lambda))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * lambda * sigma2 * x[inside])
  out
}

mp_positive_mass <- function(lambda) {
  ifelse(lambda > 1, 1 / lambda, 1)
}

shifted_mp_density <- function(x, delta, lambda, sigma2) {
  mp_density(x - delta, lambda, sigma2) / mp_positive_mass(lambda)
}

weighted_mean <- function(x, w) {
  w[!is.finite(w)] <- 0
  if (sum(w) <= 1e-10) {
    return(mean(x))
  }
  sum(w * x) / sum(w)
}

weighted_var <- function(x, w, mu) {
  w[!is.finite(w)] <- 0
  if (sum(w) <= 1e-10 || !is.finite(mu)) {
    return(var(x))
  }
  sum(w * (x - mu)^2) / sum(w)
}

moment_start <- function(vals, weights) {
  mu <- weighted_mean(vals, weights)
  v <- weighted_var(vals, weights, mu)
  lambda <- v / max(mu^2, .Machine$double.eps)
  if (!is.finite(lambda)) {
    lambda <- 0.5
  }
  active_vals <- vals[weights > median(weights)]
  if (length(active_vals) < 10) {
    active_vals <- vals
  }
  delta <- max(0, as.numeric(quantile(active_vals, 0.01, names = FALSE)) * 0.8)
  c(
    delta = delta,
    lambda = min(max(lambda, 1e-4), 0.95),
    sigma2 = max(mu - delta, .Machine$double.eps)
  )
}

fit_shifted_mp_weighted <- function(vals, weights, start) {
  weights[!is.finite(weights)] <- 0
  if (sum(weights) <= 1e-10) {
    weights <- rep(1, length(vals))
  }

  start <- c(
    delta = as.numeric(start["delta"]),
    lambda = as.numeric(start["lambda"]),
    sigma2 = as.numeric(start["sigma2"])
  )
  if (any(!is.finite(start)) || start["delta"] < 0 || start["lambda"] <= 0 || start["sigma2"] <= 0) {
    start <- moment_start(vals, weights)
  }

  active_vals <- vals[weights > median(weights)]
  if (length(active_vals) < 10) {
    active_vals <- vals
  }
  delta_upper <- max(0, as.numeric(quantile(active_vals, 0.20, names = FALSE)) * 0.999)
  mu <- weighted_mean(vals, weights)

  starts <- rbind(
    start,
    c(delta = 0, lambda = start["lambda"], sigma2 = max(mu, .Machine$double.eps)),
    c(delta = delta_upper * 0.25, lambda = start["lambda"], sigma2 = max(mu - delta_upper * 0.25, .Machine$double.eps)),
    c(delta = delta_upper * 0.75, lambda = start["lambda"], sigma2 = max(mu - delta_upper * 0.75, .Machine$double.eps)),
    c(delta = start["delta"], lambda = start["lambda"] * 0.5, sigma2 = start["sigma2"]),
    c(delta = start["delta"], lambda = min(start["lambda"] * 2, 0.95), sigma2 = start["sigma2"]),
    c(delta = start["delta"], lambda = 0.3, sigma2 = start["sigma2"]),
    c(delta = start["delta"], lambda = 0.8, sigma2 = start["sigma2"])
  )
  starts[, "delta"] <- pmin(pmax(starts[, "delta"], 0), delta_upper)
  starts[, "lambda"] <- pmin(pmax(starts[, "lambda"], 1e-4), 1)
  starts[, "sigma2"] <- pmax(starts[, "sigma2"], .Machine$double.eps)

  objective <- function(par) {
    delta <- par[1]
    lambda <- exp(par[2])
    sigma2 <- exp(par[3])
    d <- pmax(shifted_mp_density(vals, delta, lambda, sigma2), density_floor)
    -sum(weights * log(d))
  }

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = c(starts[i, "delta"], log(starts[i, c("lambda", "sigma2")])),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(0, log(1e-4), log(.Machine$double.eps)),
      upper = c(delta_upper, log(1), log(max(vals) * 100)),
      control = list(maxit = 1500)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  out <- c(
    delta = as.numeric(best$par[1]),
    lambda = as.numeric(exp(best$par[2])),
    sigma2 = as.numeric(exp(best$par[3])),
    objective = best$value,
    convergence = best$convergence
  )
  if (any(!is.finite(out)) || out["delta"] < 0 || out["lambda"] <= 0 || out["sigma2"] <= 0) {
    out <- c(start, objective = NA_real_, convergence = NA_real_)
  }
  out
}

fit_two_shifted_mp_em <- function(vals, init_resp, init_params = NULL) {
  vals <- vals[is.finite(vals) & vals > 0]
  resp <- init_resp
  resp <- resp / rowSums(resp)

  if (is.null(init_params)) {
    starts <- list(
      moment_start(vals, resp[, 1]),
      moment_start(vals, resp[, 2])
    )
  } else {
    init_params <- init_params[order(init_params$center), ]
    starts <- lapply(seq_len(2), function(i) {
      c(
        delta = init_params$delta[i],
        lambda = init_params$lambda[i],
        sigma2 = init_params$sigma2[i]
      )
    })
  }

  loglik_trace <- numeric(max_iter)
  param_trace <- matrix(NA_real_, nrow = max_iter, ncol = 8)
  colnames(param_trace) <- c("pi1", "pi2", "delta1", "delta2", "lambda1", "lambda2", "sigma1", "sigma2")

  for (iter in seq_len(max_iter)) {
    fit1 <- fit_shifted_mp_weighted(vals, resp[, 1], starts[[1]])
    fit2 <- fit_shifted_mp_weighted(vals, resp[, 2], starts[[2]])
    fits <- list(fit1, fit2)

    centers <- vapply(fits, function(fit) fit["delta"] + fit["sigma2"], numeric(1))
    ord <- order(centers)
    fits <- fits[ord]
    resp <- resp[, ord, drop = FALSE]
    starts <- lapply(fits, function(fit) c(delta = fit["delta"], lambda = fit["lambda"], sigma2 = fit["sigma2"]))

    pi_m_step <- colMeans(resp)
    d1 <- pmax(shifted_mp_density(vals, fits[[1]]["delta"], fits[[1]]["lambda"], fits[[1]]["sigma2"]), density_floor)
    d2 <- pmax(shifted_mp_density(vals, fits[[2]]["delta"], fits[[2]]["lambda"], fits[[2]]["sigma2"]), density_floor)
    weighted_density <- cbind(pi_m_step[1] * d1, pi_m_step[2] * d2)
    total_density <- pmax(rowSums(weighted_density), density_floor)
    next_resp <- weighted_density / total_density
    next_resp[!is.finite(next_resp)] <- 0
    empty <- rowSums(next_resp) <= 0
    next_resp[empty, ] <- matrix(rep(c(0.5, 0.5), sum(empty)), ncol = 2, byrow = TRUE)
    next_resp <- next_resp / rowSums(next_resp)

    loglik_trace[iter] <- sum(log(total_density))
    param_trace[iter, ] <- c(
      pi_m_step,
      fits[[1]]["delta"],
      fits[[2]]["delta"],
      fits[[1]]["lambda"],
      fits[[2]]["lambda"],
      fits[[1]]["sigma2"],
      fits[[2]]["sigma2"]
    )
    resp <- next_resp

    if (
      iter > 1 &&
        is.finite(loglik_trace[iter]) &&
        abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))
    ) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      param_trace <- param_trace[seq_len(iter), , drop = FALSE]
      break
    }
  }

  final <- param_trace[nrow(param_trace), ]
  weight <- c(final["pi1"], final["pi2"])
  delta <- c(final["delta1"], final["delta2"])
  lambda <- c(final["lambda1"], final["lambda2"])
  sigma2 <- c(final["sigma1"], final["sigma2"])

  list(
    resp = resp,
    parameters = data.frame(
      component = c("shifted_mp_low", "shifted_mp_high"),
      weight = weight,
      p_hat = weight * length(vals),
      delta = delta,
      lambda = lambda,
      sigma2 = sigma2,
      center = delta + sigma2,
      mp_lower = delta + sigma2 * (1 - sqrt(lambda))^2,
      mp_upper = delta + sigma2 * (1 + sqrt(lambda))^2,
      n_eigenvalues = length(vals),
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace)
  )
}

make_diagnostics <- function(vals, fit) {
  params <- fit$parameters
  x_grid <- seq(max(min(vals), 1e-12), max(vals, params$mp_upper, na.rm = TRUE), length.out = 25000)

  curve_df <- rbind(
    data.frame(component = params$component[1], eigenvalue = x_grid, density = params$weight[1] * shifted_mp_density(x_grid, params$delta[1], params$lambda[1], params$sigma2[1])),
    data.frame(component = params$component[2], eigenvalue = x_grid, density = params$weight[2] * shifted_mp_density(x_grid, params$delta[2], params$lambda[2], params$sigma2[2]))
  )
  mix_df <- aggregate(density ~ eigenvalue, curve_df, sum)
  dx <- c(diff(mix_df$eigenvalue), tail(diff(mix_df$eigenvalue), 1))
  mix_df$cdf <- cumsum(mix_df$density * dx)
  mix_df$cdf <- mix_df$cdf / max(mix_df$cdf)

  empirical_df <- data.frame(
    eigenvalue = sort(vals),
    empirical_cdf = seq_along(vals) / length(vals)
  )
  empirical_df$model_cdf <- approx(mix_df$eigenvalue, mix_df$cdf, xout = empirical_df$eigenvalue, rule = 2)$y

  qq_probs <- seq(0.0025, 0.9975, length.out = 600)
  qq_df <- data.frame(
    prob = qq_probs,
    empirical = as.numeric(quantile(vals, qq_probs, names = FALSE)),
    model = approx(mix_df$cdf, mix_df$eigenvalue, xout = qq_probs, ties = "ordered", rule = 2)$y
  )

  log_grid <- seq(log10(min(vals)), log10(max(vals)), length.out = 25000)
  raw_from_log <- 10^log_grid
  log_curve_df <- rbind(
    data.frame(component = params$component[1], log10_eigenvalue = log_grid, density = params$weight[1] * shifted_mp_density(raw_from_log, params$delta[1], params$lambda[1], params$sigma2[1]) * raw_from_log * log(10)),
    data.frame(component = params$component[2], log10_eigenvalue = log_grid, density = params$weight[2] * shifted_mp_density(raw_from_log, params$delta[2], params$lambda[2], params$sigma2[2]) * raw_from_log * log(10))
  )
  log_mix_df <- aggregate(density ~ log10_eigenvalue, log_curve_df, sum)

  write.csv(mix_df, file.path(out_dir, "shifted_mp_mixture_model_cdf.csv"), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, "shifted_mp_mixture_empirical_vs_model_cdf.csv"), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, "shifted_mp_mixture_qq.csv"), row.names = FALSE)

  theme_set(theme_classic(base_size = 13))

  raw_density_plot <- ggplot(data.frame(eigenvalue = vals), aes(eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = mix_df, aes(eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "Shifted-MP mixture EM initialized from GMM",
      subtitle = "Outlier component removed before fitting",
      x = "Covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  log_density_plot <- ggplot(data.frame(log10_eigenvalue = log10(vals)), aes(log10_eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = log_curve_df, aes(log10_eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = log_mix_df, aes(log10_eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "Shifted-MP mixture EM initialized from GMM, log10 view",
      subtitle = "Densities transformed with the log10 Jacobian",
      x = "log10 covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    labs(
      title = "Shifted-MP mixture CDF",
      subtitle = "Grey: empirical CDF; blue: fitted mixture CDF",
      x = "Covariance eigenvalue",
      y = "CDF"
    )

  qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.35, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = qq_lim, ylim = qq_lim) +
    labs(
      title = "Shifted-MP mixture Q-Q",
      x = "Fitted shifted-MP mixture quantile",
      y = "Empirical quantile"
    )

  ggsave(file.path(out_dir, "shifted_mp_mixture_density_raw.png"), raw_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_mixture_density_log10.png"), log_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_mixture_cdf.png"), cdf_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_mixture_qq.png"), qq_plot, width = 6, height = 6, dpi = 220)
}

gmm_assign <- read.csv(file.path(gmm_dir, "log10_gmm_assignments.csv"))
gmm_params <- read.csv(file.path(gmm_dir, "log10_parameters.csv"))
init_params_file <- file.path(init_dir, "shifted_mp_gmm_init_parameters.csv")
if (!file.exists(init_params_file)) {
  stop("Expected shifted-MP initialization parameters. Run plot_shifted_mp_gmm_initialization.R first.")
}
shifted_mp_init_params <- read.csv(init_params_file)

if (!all(c("log_g1", "log_g2", "log_outlier") %in% names(gmm_assign))) {
  stop("Expected log_g1, log_g2, and log_outlier columns. Re-run fit_gmm_spectrum.R first.")
}

keep <- gmm_assign$log_outlier < 0.5
bulk <- gmm_assign[keep, ]
vals <- bulk$eigenvalue
init_resp <- as.matrix(bulk[, c("log_g1", "log_g2")])
init_resp <- init_resp / rowSums(init_resp)

fit <- fit_two_shifted_mp_em(vals, init_resp, shifted_mp_init_params)

assign_df <- data.frame(
  n_cells = bulk$n_cells,
  rank = bulk$rank,
  eigenvalue = vals,
  log10_eigenvalue = log10(vals),
  gmm_component = bulk$log_component,
  shifted_mp_low = fit$resp[, 1],
  shifted_mp_high = fit$resp[, 2]
)
assign_df$mp_component <- fit$parameters$component[max.col(fit$resp, ties.method = "first")]

initial_summary <- data.frame(
  component = c("gmm_g1", "gmm_g2", "removed_outlier"),
  p_hat = c(sum(init_resp[, 1]), sum(init_resp[, 2]), sum(!keep)),
  source = "log10_gmm_plus_outlier"
)

write.csv(gmm_params, file.path(out_dir, "gmm_initial_parameters.csv"), row.names = FALSE)
write.csv(shifted_mp_init_params, file.path(out_dir, "shifted_mp_starting_parameters.csv"), row.names = FALSE)
write.csv(initial_summary, file.path(out_dir, "gmm_initial_bulk_counts.csv"), row.names = FALSE)
write.csv(fit$parameters, file.path(out_dir, "shifted_mp_mixture_parameters.csv"), row.names = FALSE)
write.csv(assign_df, file.path(out_dir, "shifted_mp_mixture_assignments.csv"), row.names = FALSE)
write.csv(fit$loglik_trace, file.path(out_dir, "shifted_mp_mixture_loglik.csv"), row.names = FALSE)

make_diagnostics(vals, fit)

print(initial_summary)
print(fit$parameters)
message("Wrote shifted-MP mixture EM outputs to ", out_dir)
