library(ggplot2)

out_dir <- "outputs/exploration/test_single_shifted_mp_fit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

n <- 1200
p <- 360
lambda_true <- p / n
sigma2_true <- 0.18
delta_true <- 0.04

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

simulate_shifted_gram_eigenvalues <- function(n, p, sigma2, delta) {
  x <- matrix(rnorm(p * n, sd = sqrt(sigma2)), nrow = p, ncol = n)
  gram <- tcrossprod(x) / n
  shifted_gram <- gram + diag(delta, nrow = p)
  sort(eigen(shifted_gram, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
}

simulate_shifted_gram_eigenvalues_fast <- function(n, p, sigma2, delta) {
  x <- matrix(rnorm(p * n, sd = sqrt(sigma2)), nrow = p, ncol = n)
  gram <- tcrossprod(x) / n
  sort(eigen(gram, symmetric = TRUE, only.values = TRUE)$values + delta, decreasing = TRUE)
}

fit_shifted_mp <- function(vals) {
  vals <- vals[is.finite(vals) & vals > 0]

  objective <- function(par) {
    delta <- par[1]
    lambda <- exp(par[2])
    sigma2 <- exp(par[3])
    d <- pmax(shifted_mp_density(vals, delta, lambda, sigma2), density_floor)
    -sum(log(d))
  }

  q <- quantile(vals, c(0.01, 0.05, 0.50, 0.95), names = FALSE)
  starts <- expand.grid(
    delta = c(0, q[1] * 0.5, q[1] * 0.9, q[2] * 0.8),
    lambda = c(0.05, 0.15, 0.3, 0.6, 0.9),
    sigma2 = c(mean(vals), median(vals), q[3], diff(range(vals)) / 4)
  )
  starts$delta <- pmin(pmax(starts$delta, 0), min(vals) * 0.999)
  starts$sigma2 <- pmax(starts$sigma2 - starts$delta, .Machine$double.eps)

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = c(starts$delta[i], log(starts$lambda[i]), log(starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(0, log(1e-4), log(.Machine$double.eps)),
      upper = c(min(vals) * 0.999, log(1), log(max(vals) * 100)),
      control = list(maxit = 2000)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  delta <- best$par[1]
  lambda <- exp(best$par[2])
  sigma2 <- exp(best$par[3])
  data.frame(
    delta = delta,
    lambda = lambda,
    sigma2 = sigma2,
    mp_lower = delta + sigma2 * (1 - sqrt(lambda))^2,
    mp_upper = delta + sigma2 * (1 + sqrt(lambda))^2,
    neg_loglik = best$value,
    convergence = best$convergence
  )
}

make_diagnostics <- function(vals, fit) {
  x_grid <- seq(min(vals), max(vals), length.out = 20000)
  density_df <- data.frame(
    eigenvalue = x_grid,
    density = shifted_mp_density(x_grid, fit$delta, fit$lambda, fit$sigma2)
  )
  dx <- c(diff(x_grid), tail(diff(x_grid), 1))
  density_df$cdf <- cumsum(density_df$density * dx)
  density_df$cdf <- density_df$cdf / max(density_df$cdf)

  empirical_df <- data.frame(
    eigenvalue = sort(vals),
    empirical_cdf = seq_along(vals) / length(vals)
  )
  empirical_df$model_cdf <- approx(density_df$eigenvalue, density_df$cdf, xout = empirical_df$eigenvalue, rule = 2)$y

  probs <- seq(0.0025, 0.9975, length.out = 600)
  qq_df <- data.frame(
    prob = probs,
    empirical = as.numeric(quantile(vals, probs, names = FALSE)),
    model = approx(density_df$cdf, density_df$eigenvalue, xout = probs, ties = "ordered", rule = 2)$y
  )

  density_plot <- ggplot(data.frame(eigenvalue = vals), aes(eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 90, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = density_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 1) +
    geom_vline(xintercept = c(fit$mp_lower, fit$mp_upper), color = "#C43C39", linetype = "dashed") +
    labs(
      title = "Single shifted-MP fit on simulated shifted Gram spectrum",
      x = "Eigenvalue",
      y = "Density"
    )

  cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    labs(
      title = "Single shifted-MP CDF check",
      subtitle = "Grey: empirical CDF; blue: fitted shifted-MP CDF",
      x = "Eigenvalue",
      y = "CDF"
    )

  qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.35, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = qq_lim, ylim = qq_lim) +
    labs(
      title = "Single shifted-MP Q-Q check",
      x = "Fitted shifted-MP quantile",
      y = "Empirical quantile"
    )

  write.csv(density_df, file.path(out_dir, "single_shifted_mp_model_cdf.csv"), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, "single_shifted_mp_empirical_vs_model_cdf.csv"), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, "single_shifted_mp_qq.csv"), row.names = FALSE)

  ggsave(file.path(out_dir, "single_shifted_mp_density.png"), density_plot, width = 8, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "single_shifted_mp_cdf.png"), cdf_plot, width = 8, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "single_shifted_mp_qq.png"), qq_plot, width = 6, height = 6, dpi = 220)
}

one_fit_row <- function(scenario_id, n, p, sigma2, delta) {
  vals <- simulate_shifted_gram_eigenvalues_fast(n, p, sigma2, delta)
  fit <- fit_shifted_mp(vals)
  lambda <- p / n
  data.frame(
    scenario_id = scenario_id,
    n = n,
    p = p,
    true_delta = delta,
    fit_delta = fit$delta,
    true_lambda = lambda,
    fit_lambda = fit$lambda,
    true_sigma2 = sigma2,
    fit_sigma2 = fit$sigma2,
    true_lower = delta + sigma2 * (1 - sqrt(lambda))^2,
    fit_lower = fit$mp_lower,
    true_upper = delta + sigma2 * (1 + sqrt(lambda))^2,
    fit_upper = fit$mp_upper,
    neg_loglik = fit$neg_loglik,
    convergence = fit$convergence
  )
}

manual_shift_row <- function(shift_id, base_vals, lambda, sigma2, delta) {
  vals <- base_vals + delta
  fit <- fit_shifted_mp(vals)
  data.frame(
    shift_id = shift_id,
    true_delta = delta,
    fit_delta = fit$delta,
    true_lambda = lambda,
    fit_lambda = fit$lambda,
    true_sigma2 = sigma2,
    fit_sigma2 = fit$sigma2,
    true_lower = delta + sigma2 * (1 - sqrt(lambda))^2,
    fit_lower = fit$mp_lower,
    true_upper = delta + sigma2 * (1 + sqrt(lambda))^2,
    fit_upper = fit$mp_upper,
    neg_loglik = fit$neg_loglik,
    convergence = fit$convergence
  )
}

eigvals <- simulate_shifted_gram_eigenvalues(n, p, sigma2_true, delta_true)
fit <- fit_shifted_mp(eigvals)

truth <- data.frame(
  parameter = c("delta", "lambda", "sigma2", "mp_lower", "mp_upper"),
  true = c(
    delta_true,
    lambda_true,
    sigma2_true,
    delta_true + sigma2_true * (1 - sqrt(lambda_true))^2,
    delta_true + sigma2_true * (1 + sqrt(lambda_true))^2
  ),
  fitted = c(fit$delta, fit$lambda, fit$sigma2, fit$mp_lower, fit$mp_upper)
)

write.csv(data.frame(eigenvalue = eigvals), file.path(out_dir, "simulated_shifted_gram_eigenvalues.csv"), row.names = FALSE)
write.csv(fit, file.path(out_dir, "single_shifted_mp_fit.csv"), row.names = FALSE)
write.csv(truth, file.path(out_dir, "single_shifted_mp_truth_vs_fit.csv"), row.names = FALSE)

scenarios <- data.frame(
  scenario_id = seq_len(4),
  n = c(800, 800, 800, 700),
  p = c(240, 240, 400, 350),
  sigma2 = c(0.18, 0.18, 0.10, 0.25),
  delta = c(0.00, 0.04, 0.04, 0.08)
)
robustness <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  with(scenarios[i, ], one_fit_row(scenario_id, n, p, sigma2, delta))
}))
write.csv(robustness, file.path(out_dir, "single_shifted_mp_robustness_sweep.csv"), row.names = FALSE)

manual_n <- 1200
manual_p <- 360
manual_lambda <- manual_p / manual_n
manual_sigma2 <- 0.18
base_vals <- simulate_shifted_gram_eigenvalues_fast(manual_n, manual_p, manual_sigma2, delta = 0)
manual_shifts <- c(0, 0.02, 0.05, 0.10, 0.20, 0.40)
manual_shift_sweep <- do.call(rbind, lapply(seq_along(manual_shifts), function(i) {
  manual_shift_row(i, base_vals, manual_lambda, manual_sigma2, manual_shifts[i])
}))
write.csv(manual_shift_sweep, file.path(out_dir, "single_shifted_mp_manual_shift_sweep.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))
make_diagnostics(eigvals, fit)

print(truth)
print(robustness)
print(manual_shift_sweep)
message("Wrote single shifted-MP simulation fit outputs to ", out_dir)
