library(ggplot2)

out_dir <- "outputs/exploration/eigenvalue_bulk_shape"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  jurkat_hvg = list(
    label = "Jurkat HVG",
    noise_file = "outputs/exploration/jurkat_hvg_noise_mp_by_cells/jurkat_hvg_noise_eigenvalues_by_cells.csv",
    fit_file = "outputs/exploration/jurkat_hvg_noise_mp_by_cells/jurkat_hvg_noise_threshold_mp_fit_by_cells.csv"
  ),
  k562_hvg = list(
    label = "K562 HVG",
    noise_file = "outputs/exploration/k562_320k_hvg_noise_mp_by_cells/k562_hvg_noise_eigenvalues_by_cells.csv",
    fit_file = "outputs/exploration/k562_320k_hvg_noise_mp_by_cells/k562_hvg_noise_threshold_mp_fit_by_cells.csv"
  )
)

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_quantiles <- function(gamma, sigma2, probs, grid_n = 10000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 1e-12), upper - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_gamma_mle <- function(x) {
  objective <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    -sum(dgamma(x, shape = shape, scale = scale, log = TRUE))
  }
  start <- c(log(mean(x)^2 / var(x)), log(var(x) / mean(x)))
  fit <- optim(start, objective, method = "Nelder-Mead", control = list(maxit = 1000))
  c(shape = exp(fit$par[1]), scale = exp(fit$par[2]))
}

fit_lognormal_mixture <- function(x, n_iter = 300) {
  z <- log(x)
  qs <- quantile(z, c(0.35, 0.75), names = FALSE)
  mu <- c(qs[1], qs[2])
  sigma <- rep(sd(z), 2)
  weight <- c(0.65, 0.35)

  for (iter in seq_len(n_iter)) {
    d1 <- weight[1] * dnorm(z, mu[1], sigma[1])
    d2 <- weight[2] * dnorm(z, mu[2], sigma[2])
    total <- pmax(d1 + d2, .Machine$double.xmin)
    r1 <- d1 / total
    r2 <- 1 - r1
    weight <- c(mean(r1), mean(r2))
    mu <- c(sum(r1 * z) / sum(r1), sum(r2 * z) / sum(r2))
    sigma <- c(
      sqrt(sum(r1 * (z - mu[1])^2) / sum(r1)),
      sqrt(sum(r2 * (z - mu[2])^2) / sum(r2))
    )
    sigma <- pmax(sigma, 1e-4)
  }

  ord <- order(mu)
  list(weight = weight[ord], mu = mu[ord], sigma = sigma[ord])
}

mixture_quantiles <- function(fit, probs, x) {
  grid <- exp(seq(log(min(x) * 0.25), log(max(x) * 4), length.out = 30000))
  cdf <- fit$weight[1] * plnorm(grid, fit$mu[1], fit$sigma[1]) +
    fit$weight[2] * plnorm(grid, fit$mu[2], fit$sigma[2])
  approx(cdf, grid, xout = probs, ties = "ordered", rule = 2)$y
}

quantile_rmse <- function(empirical_q, fitted_q) {
  sqrt(mean((log(empirical_q) - log(fitted_q))^2))
}

diagnose_dataset <- function(dataset_name, spec) {
  noise_df <- read.csv(spec$noise_file)
  fit_df <- read.csv(spec$fit_file)
  n_cells <- max(noise_df$n_cells)

  x_all <- noise_df$eigenvalue[noise_df$n_cells == n_cells]
  x <- sort(x_all[x_all > 1e-10])
  fit_row <- fit_df[fit_df$n_cells == n_cells, ]
  probs <- seq(0.01, 0.99, by = 0.01)
  empirical_q <- as.numeric(quantile(x, probs, names = FALSE))

  gamma_fit <- fit_gamma_mle(x)
  lognormal_mean <- mean(log(x))
  lognormal_sd <- sd(log(x))
  mixture_fit <- fit_lognormal_mixture(x)

  model_q <- data.frame(
    prob = rep(probs, 4),
    empirical = rep(empirical_q, 4),
    fitted = c(
      mp_quantiles(fit_row$mp_gamma, fit_row$mp_sigma2, probs),
      qgamma(probs, shape = gamma_fit["shape"], scale = gamma_fit["scale"]),
      qlnorm(probs, meanlog = lognormal_mean, sdlog = lognormal_sd),
      mixture_quantiles(mixture_fit, probs, x)
    ),
    model = rep(
      c("Marchenko-Pastur", "Gamma", "Lognormal", "2-lognormal mixture"),
      each = length(probs)
    ),
    dataset = spec$label,
    n_cells = n_cells
  )
  model_q$log_residual <- log(model_q$empirical) - log(model_q$fitted)

  metrics <- aggregate(
    cbind(log_residual) ~ dataset + n_cells + model,
    model_q,
    function(v) sqrt(mean(v^2))
  )
  names(metrics)[names(metrics) == "log_residual"] <- "log_quantile_rmse"

  tail_probs <- c(0.5, 0.9, 0.95, 0.99)
  empirical_tail <- as.numeric(quantile(x, tail_probs, names = FALSE))
  mp_tail <- mp_quantiles(fit_row$mp_gamma, fit_row$mp_sigma2, tail_probs)
  tail_df <- data.frame(
    dataset = spec$label,
    n_cells = n_cells,
    prob = tail_probs,
    empirical = empirical_tail,
    mp = mp_tail,
    empirical_over_mp = empirical_tail / mp_tail
  )

  summary <- data.frame(
    dataset = spec$label,
    n_cells = n_cells,
    n_positive = length(x),
    n_zero = sum(x_all <= 1e-10),
    min_positive = min(x),
    median = median(x),
    mean = mean(x),
    q90 = quantile(x, 0.90, names = FALSE),
    q99 = quantile(x, 0.99, names = FALSE),
    max = max(x),
    mp_gamma = fit_row$mp_gamma,
    mp_sigma2 = fit_row$mp_sigma2,
    mp_upper = fit_row$mp_upper,
    empirical_q99_over_mp_q99 = tail_df$empirical_over_mp[tail_df$prob == 0.99],
    mixture_upper_weight = mixture_fit$weight[2],
    mixture_lower_median = exp(mixture_fit$mu[1]),
    mixture_upper_median = exp(mixture_fit$mu[2])
  )

  list(quantiles = model_q, metrics = metrics, tail = tail_df, summary = summary)
}

results <- Map(diagnose_dataset, names(datasets), datasets)

quantile_df <- do.call(rbind, lapply(results, `[[`, "quantiles"))
metrics_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
tail_df <- do.call(rbind, lapply(results, `[[`, "tail"))
summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))

write.csv(metrics_df, file.path(out_dir, "bulk_model_quantile_rmse.csv"), row.names = FALSE)
write.csv(tail_df, file.path(out_dir, "mp_tail_quantile_ratios.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "bulk_shape_summary.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

qq_plot <- ggplot(quantile_df, aes(fitted, empirical, color = model)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ dataset, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Noise-bulk eigenvalue quantile diagnostics",
    subtitle = "Full-cell positive noise eigenvalues after the scree-elbow spike cut",
    x = "Model quantile",
    y = "Empirical quantile",
    color = "Model"
  )

residual_plot <- ggplot(quantile_df, aes(prob, log_residual, color = model)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ dataset) +
  labs(
    title = "Log-quantile residuals by model",
    x = "Quantile probability",
    y = "log(empirical / fitted)",
    color = "Model"
  )

ggsave(
  file.path(out_dir, "bulk_model_qq.png"),
  qq_plot,
  width = 9,
  height = 5.5,
  dpi = 220
)
ggsave(
  file.path(out_dir, "bulk_model_log_quantile_residuals.png"),
  residual_plot,
  width = 9,
  height = 5.5,
  dpi = 220
)

print(summary_df)
print(metrics_df[order(metrics_df$dataset, metrics_df$log_quantile_rmse), ])
message("Wrote eigenvalue bulk diagnostics to ", out_dir)
