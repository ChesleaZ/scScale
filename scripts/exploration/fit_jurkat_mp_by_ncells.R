library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  in_file <- args[[1]]
  out_dir <- args[[2]]
  dataset_label <- args[[3]]
} else {
  in_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
  out_dir <- "outputs/exploration/jurkat_mp_by_ncells"
  dataset_label <- "Jurkat"
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_probs <- seq(0.05, 0.95, by = 0.025)
grid_n <- 20000

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_positive_quantiles <- function(gamma, sigma2, probs, grid_n = 20000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 0) + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_mp_fixed_gamma <- function(vals, gamma) {
  positive_vals <- vals[is.finite(vals) & vals > 1e-10]
  empirical_q <- as.numeric(quantile(positive_vals, fit_probs, names = FALSE))

  objective <- function(log_sigma2) {
    sigma2 <- exp(log_sigma2)
    fitted_q <- mp_positive_quantiles(gamma, sigma2, fit_probs, grid_n = grid_n)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- log(c(median(positive_vals), mean(positive_vals), quantile(positive_vals, 0.75)))
  fits <- lapply(starts, function(start) {
    optim(start, objective, method = "Brent", lower = log(min(positive_vals) / 100), upper = log(max(positive_vals) * 100))
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  sigma2 <- unname(exp(best$par))
  c(
    gamma = gamma,
    sigma2 = sigma2,
    mp_lower = sigma2 * (1 - sqrt(gamma))^2,
    mp_upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = best$value,
    convergence = best$convergence
  )
}

fit_mp_free_gamma <- function(vals, gamma_start) {
  positive_vals <- vals[is.finite(vals) & vals > 1e-10]
  empirical_q <- as.numeric(quantile(positive_vals, fit_probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_positive_quantiles(gamma, sigma2, fit_probs, grid_n = grid_n)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = unique(pmax(c(gamma_start, 0.25, 0.5, 1, 2, 4), 1e-4)),
    sigma2 = c(median(positive_vals), mean(positive_vals), quantile(positive_vals, 0.75))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = log(c(1e-4, min(positive_vals) / 100)),
      upper = log(c(25, max(positive_vals) * 100)),
      control = list(maxit = 1200)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma <- unname(exp(best$par[1]))
  sigma2 <- unname(exp(best$par[2]))
  c(
    gamma = gamma,
    sigma2 = sigma2,
    mp_lower = sigma2 * (1 - sqrt(gamma))^2,
    mp_upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = best$value,
    convergence = best$convergence
  )
}

make_curve <- function(n_cells, vals, fit, fit_type) {
  positive_vals <- vals[is.finite(vals) & vals > 1e-10]
  x <- seq(min(positive_vals), quantile(positive_vals, 0.995, names = FALSE), length.out = 2000)
  data.frame(
    n_cells = n_cells,
    fit_type = fit_type,
    eigenvalue = x,
    density = mp_density(x, fit["gamma"], fit["sigma2"])
  )
}

eig_df <- read.csv(in_file)
n_values <- sort(unique(eig_df$n_cells))
p <- max(eig_df$rank)

fit_rows <- list()
curve_rows <- list()
for (n_cells in n_values) {
  vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  vals <- vals[is.finite(vals)]
  gamma_theory <- p / n_cells

  fixed_fit <- fit_mp_fixed_gamma(vals, gamma_theory)
  free_fit <- fit_mp_free_gamma(vals, gamma_theory)

  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    p = p,
    n_eigenvalues = length(vals),
    n_positive = sum(vals > 1e-10),
    n_zero = sum(vals <= 1e-10),
    zero_fraction = mean(vals <= 1e-10),
    gamma_theory = gamma_theory,
    fit_type = "fixed_gamma",
    gamma = fixed_fit["gamma"],
    sigma2 = fixed_fit["sigma2"],
    mp_lower = fixed_fit["mp_lower"],
    mp_upper = fixed_fit["mp_upper"],
    objective = fixed_fit["objective"],
    convergence = fixed_fit["convergence"]
  )
  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    p = p,
    n_eigenvalues = length(vals),
    n_positive = sum(vals > 1e-10),
    n_zero = sum(vals <= 1e-10),
    zero_fraction = mean(vals <= 1e-10),
    gamma_theory = gamma_theory,
    fit_type = "free_gamma",
    gamma = free_fit["gamma"],
    sigma2 = free_fit["sigma2"],
    mp_lower = free_fit["mp_lower"],
    mp_upper = free_fit["mp_upper"],
    objective = free_fit["objective"],
    convergence = free_fit["convergence"]
  )

  curve_rows[[length(curve_rows) + 1]] <- make_curve(n_cells, vals, fixed_fit, "fixed_gamma")
  curve_rows[[length(curve_rows) + 1]] <- make_curve(n_cells, vals, free_fit, "free_gamma")
}

fit_df <- do.call(rbind, fit_rows)
curve_df <- do.call(rbind, curve_rows)

write.csv(fit_df, file.path(out_dir, "jurkat_mp_by_ncells_parameters.csv"), row.names = FALSE)
write.csv(curve_df, file.path(out_dir, "jurkat_mp_by_ncells_curves.csv"), row.names = FALSE)

free_fit_df <- fit_df[fit_df$fit_type == "free_gamma", ]
free_curve_df <- curve_df[curve_df$fit_type == "free_gamma", ]

plot_df <- fit_df
plot_df$fit_type <- factor(plot_df$fit_type, levels = c("fixed_gamma", "free_gamma"))

param_long <- rbind(
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "gamma", value = plot_df$gamma),
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "sigma2", value = plot_df$sigma2),
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "mp_lower", value = plot_df$mp_lower),
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "mp_upper", value = plot_df$mp_upper),
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "objective", value = plot_df$objective),
  data.frame(n_cells = plot_df$n_cells, fit_type = plot_df$fit_type, parameter = "zero_fraction", value = plot_df$zero_fraction)
)

theme_set(theme_classic(base_size = 13))

param_plot <- ggplot(param_long, aes(n_cells, value, color = fit_type)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_line(
    data = data.frame(n_cells = n_values, parameter = "gamma", value = p / n_values),
    aes(n_cells, value),
    inherit.aes = FALSE,
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_log10(breaks = n_values, labels = format(n_values, big.mark = ",")) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_label, " one-mode MP fit across cell counts"),
    subtitle = "Dashed grey line in gamma panel is theoretical p / n",
    x = "Number of cells",
    y = "Parameter value",
    color = "Fit"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "jurkat_mp_parameters_by_ncells.png"), param_plot, width = 10, height = 8, dpi = 220)

positive_df <- eig_df[eig_df$eigenvalue > 1e-10, ]
positive_df$n_cells_label <- factor(
  positive_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)
curve_df$n_cells_label <- factor(
  curve_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)

density_plot <- ggplot(positive_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 90, fill = "grey80", color = "white", linewidth = 0.1) +
  geom_line(data = curve_df, aes(eigenvalue, density, color = fit_type), linewidth = 0.8) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    title = paste0(dataset_label, " positive eigenvalues with one-mode MP fits"),
    subtitle = "Positive continuous bulk only; robust quantile fit from 5th to 95th percentiles",
    x = "Covariance eigenvalue",
    y = "Density",
    color = "Fit"
  )

ggsave(file.path(out_dir, "jurkat_mp_density_by_ncells.png"), density_plot, width = 10, height = 8, dpi = 220)

free_density_plot <- ggplot(positive_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 90, fill = "grey80", color = "white", linewidth = 0.1) +
  geom_line(data = free_curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.9) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    title = paste0(dataset_label, " positive eigenvalues with free-gamma MP fits"),
    subtitle = "Positive continuous bulk only; robust quantile fit from 5th to 95th percentiles",
    x = "Covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "jurkat_free_gamma_mp_density_by_ncells.png"), free_density_plot, width = 10, height = 8, dpi = 220)

make_free_diagnostics <- function(n_cells, vals, fit) {
  positive_vals <- sort(vals[is.finite(vals) & vals > 1e-10])
  x_grid <- seq(min(positive_vals), max(positive_vals), length.out = 40000)
  d <- mp_density(x_grid, fit$gamma, fit$sigma2)
  dx <- c(diff(x_grid), tail(diff(x_grid), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)

  cdf_df <- data.frame(
    n_cells = n_cells,
    eigenvalue = positive_vals,
    empirical_cdf = seq_along(positive_vals) / length(positive_vals),
    model_cdf = approx(x_grid, cdf, xout = positive_vals, rule = 2)$y
  )

  probs <- seq(0.005, 0.995, length.out = 600)
  qq_df <- data.frame(
    n_cells = n_cells,
    prob = probs,
    empirical = as.numeric(quantile(positive_vals, probs, names = FALSE)),
    model = approx(cdf, x_grid, xout = probs, ties = "ordered", rule = 2)$y
  )
  list(cdf = cdf_df, qq = qq_df)
}

diag_list <- lapply(n_values, function(n_cells) {
  vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  fit <- free_fit_df[free_fit_df$n_cells == n_cells, ]
  make_free_diagnostics(n_cells, vals, fit)
})
cdf_df <- do.call(rbind, lapply(diag_list, `[[`, "cdf"))
qq_df <- do.call(rbind, lapply(diag_list, `[[`, "qq"))

cdf_df$n_cells_label <- factor(
  cdf_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)
qq_df$n_cells_label <- factor(
  qq_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)

write.csv(cdf_df, file.path(out_dir, "jurkat_free_gamma_mp_empirical_vs_model_cdf.csv"), row.names = FALSE)
write.csv(qq_df, file.path(out_dir, "jurkat_free_gamma_mp_qq.csv"), row.names = FALSE)

free_cdf_plot <- ggplot(cdf_df, aes(eigenvalue)) +
  geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.7) +
  geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.85) +
  facet_wrap(~ n_cells_label, scales = "free_x", ncol = 2) +
  labs(
    title = paste0(dataset_label, " free-gamma MP CDF diagnostics"),
    subtitle = "Grey: empirical CDF of positive eigenvalues; blue: fitted MP CDF",
    x = "Covariance eigenvalue",
    y = "CDF"
  )

ggsave(file.path(out_dir, "jurkat_free_gamma_mp_cdf_by_ncells.png"), free_cdf_plot, width = 10, height = 8, dpi = 220)

free_qq_plot <- ggplot(qq_df, aes(model, empirical)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
  geom_point(size = 1, alpha = 0.65, color = "#1F78B4") +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_label, " free-gamma MP Q-Q diagnostics"),
    x = "Fitted free-gamma MP quantile",
    y = "Empirical positive eigenvalue quantile"
  )

ggsave(file.path(out_dir, "jurkat_free_gamma_mp_qq_by_ncells.png"), free_qq_plot, width = 10, height = 8, dpi = 220)

print(fit_df)
message("Wrote MP-by-n outputs to ", out_dir)
