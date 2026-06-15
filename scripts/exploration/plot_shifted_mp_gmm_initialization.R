library(ggplot2)

gmm_dir <- "outputs/exploration/gmm_spectrum"
out_dir <- "outputs/exploration/shifted_mp_gmm_initialization"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

shifted_mp_quantiles <- function(delta, lambda, sigma2, probs, grid_n = 12000) {
  lower <- delta + sigma2 * (1 - sqrt(lambda))^2
  upper <- delta + sigma2 * (1 + sqrt(lambda))^2
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    return(rep(NA_real_, length(probs)))
  }
  x <- seq(max(lower, 1e-12), upper, length.out = grid_n)
  d <- shifted_mp_density(x, delta, lambda, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  if (!is.finite(max(cdf)) || max(cdf) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

initialize_shifted_mp_from_log_gaussian <- function(mean_log10, sd_log10, weight, component) {
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  target_q <- 10^(mean_log10 + qnorm(probs) * sd_log10)
  target_center <- 10^mean_log10

  objective <- function(par) {
    delta <- par[1]
    lambda <- exp(par[2])
    sigma2 <- exp(par[3])
    q <- shifted_mp_quantiles(delta, lambda, sigma2, probs)
    if (any(!is.finite(q)) || any(q <= 0)) {
      return(1e100)
    }
    quantile_loss <- mean((log(q) - log(target_q))^2)
    center_loss <- (log(delta + sigma2) - log(target_center))^2
    quantile_loss + 0.25 * center_loss
  }

  delta_upper <- max(0, target_q[1] * 0.999)
  starts <- expand.grid(
    delta = c(0, target_q[1] * 0.25, target_q[1] * 0.50, target_q[1] * 0.90),
    lambda = c(0.05, 0.15, 0.30, 0.60, 0.90),
    sigma2 = c(target_center, target_q[3], mean(target_q[c(2, 4)]))
  )
  starts$delta <- pmin(pmax(starts$delta, 0), delta_upper)
  starts$sigma2 <- pmax(starts$sigma2 - starts$delta, .Machine$double.eps)
  start_obj <- vapply(seq_len(nrow(starts)), function(i) {
    objective(c(starts$delta[i], log(starts$lambda[i]), log(starts$sigma2[i])))
  }, numeric(1))
  starts <- starts[is.finite(start_obj), , drop = FALSE]
  if (nrow(starts) == 0) {
    stop("No finite shifted-MP starts for component ", component)
  }

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = c(starts$delta[i], log(starts$lambda[i]), log(starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(0, log(1e-4), log(.Machine$double.eps)),
      upper = c(delta_upper, log(1), log(max(target_q) * 100)),
      control = list(maxit = 1200)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  delta <- best$par[1]
  lambda <- exp(best$par[2])
  sigma2 <- exp(best$par[3])

  data.frame(
    component = component,
    weight = weight,
    delta = delta,
    lambda = lambda,
    sigma2 = sigma2,
    center = delta + sigma2,
    mp_lower = delta + sigma2 * (1 - sqrt(lambda))^2,
    mp_upper = delta + sigma2 * (1 + sqrt(lambda))^2,
    gaussian_mean_log10 = mean_log10,
    gaussian_sd_log10 = sd_log10,
    gaussian_q05 = target_q[1],
    gaussian_q25 = target_q[2],
    gaussian_q50 = target_q[3],
    gaussian_q75 = target_q[4],
    gaussian_q95 = target_q[5],
    objective = best$value,
    convergence = best$convergence
  )
}

make_diagnostics <- function(vals, params) {
  x_grid <- seq(max(min(vals), 1e-12), max(vals, params$mp_upper, na.rm = TRUE), length.out = 25000)

  curve_df <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    data.frame(
      component = params$component[i],
      eigenvalue = x_grid,
      density = params$weight[i] * shifted_mp_density(x_grid, params$delta[i], params$lambda[i], params$sigma2[i])
    )
  }))
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
  log_curve_df <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    data.frame(
      component = params$component[i],
      log10_eigenvalue = log_grid,
      density = params$weight[i] *
        shifted_mp_density(raw_from_log, params$delta[i], params$lambda[i], params$sigma2[i]) *
        raw_from_log * log(10)
    )
  }))
  log_mix_df <- aggregate(density ~ log10_eigenvalue, log_curve_df, sum)

  write.csv(mix_df, file.path(out_dir, "shifted_mp_gmm_init_model_cdf.csv"), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, "shifted_mp_gmm_init_empirical_vs_model_cdf.csv"), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, "shifted_mp_gmm_init_qq.csv"), row.names = FALSE)

  theme_set(theme_classic(base_size = 13))

  raw_density_plot <- ggplot(data.frame(eigenvalue = vals), aes(eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = mix_df, aes(eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "Shifted-MP initialized from fitted log10-Gaussian modes",
      subtitle = "No EM; MP parameters are matched to Gaussian-implied quantiles",
      x = "Covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  log_density_plot <- ggplot(data.frame(log10_eigenvalue = log10(vals)), aes(log10_eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = log_curve_df, aes(log10_eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = log_mix_df, aes(log10_eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "Shifted-MP initialized from fitted log10-Gaussian modes, log10 view",
      subtitle = "No EM; densities transformed with the log10 Jacobian",
      x = "log10 covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    labs(
      title = "Shifted-MP GMM initialization CDF",
      subtitle = "Grey: empirical CDF; blue: initialized shifted-MP mixture CDF",
      x = "Covariance eigenvalue",
      y = "CDF"
    )

  qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.35, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = qq_lim, ylim = qq_lim) +
    labs(
      title = "Shifted-MP GMM initialization Q-Q",
      x = "Initialized shifted-MP mixture quantile",
      y = "Empirical quantile"
    )

  ggsave(file.path(out_dir, "shifted_mp_gmm_init_density_raw.png"), raw_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_gmm_init_density_log10.png"), log_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_gmm_init_cdf.png"), cdf_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "shifted_mp_gmm_init_qq.png"), qq_plot, width = 6, height = 6, dpi = 220)
}

gmm_assign <- read.csv(file.path(gmm_dir, "log10_gmm_assignments.csv"))
gmm_params <- read.csv(file.path(gmm_dir, "log10_parameters.csv"))

if (!all(c("log_g1", "log_g2", "log_outlier") %in% names(gmm_assign))) {
  stop("Expected log_g1, log_g2, and log_outlier columns. Re-run fit_gmm_spectrum.R first.")
}

mode_params <- subset(gmm_params, component %in% c("g1", "g2"))
mode_params <- mode_params[order(mode_params$mean), ]
mode_params$weight <- mode_params$weight / sum(mode_params$weight)

params <- do.call(rbind, lapply(seq_len(nrow(mode_params)), function(i) {
  initialize_shifted_mp_from_log_gaussian(
    mean_log10 = mode_params$mean[i],
    sd_log10 = mode_params$sd[i],
    weight = mode_params$weight[i],
    component = paste0("init_shifted_mp_", mode_params$component[i])
  )
}))

keep <- gmm_assign$log_outlier < 0.5
bulk <- gmm_assign[keep, ]
vals <- bulk$eigenvalue

assign_df <- data.frame(
  n_cells = bulk$n_cells,
  rank = bulk$rank,
  eigenvalue = vals,
  log10_eigenvalue = log10(vals),
  gmm_component = bulk$log_component,
  log_g1 = bulk$log_g1,
  log_g2 = bulk$log_g2
)

initial_summary <- data.frame(
  component = c("gmm_g1", "gmm_g2", "removed_outlier"),
  p_hat = c(sum(bulk$log_g1), sum(bulk$log_g2), sum(!keep)),
  source = "log10_gmm_plus_outlier"
)

write.csv(gmm_params, file.path(out_dir, "gmm_initial_parameters.csv"), row.names = FALSE)
write.csv(initial_summary, file.path(out_dir, "gmm_initial_bulk_counts.csv"), row.names = FALSE)
write.csv(params, file.path(out_dir, "shifted_mp_gmm_init_parameters.csv"), row.names = FALSE)
write.csv(assign_df, file.path(out_dir, "shifted_mp_gmm_init_assignments.csv"), row.names = FALSE)

make_diagnostics(vals, params)

print(initial_summary)
print(params)
message("Wrote shifted-MP GMM initialization outputs to ", out_dir)
