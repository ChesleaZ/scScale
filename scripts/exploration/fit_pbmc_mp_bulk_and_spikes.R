#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: fit_pbmc_mp_bulk_and_spikes.R <eigenvalues_csv> <out_dir> [dataset_name]", call. = FALSE)
}

eig_file <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

fit_probs <- seq(0.05, 0.80, by = 0.025)
grid_n <- 20000

mp_density <- function(x, gamma, sigma2) {
  if (!is.finite(gamma) || !is.finite(sigma2) || gamma <= 0 || sigma2 <= 0) {
    return(rep(0, length(x)))
  }
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_positive_mass <- function(gamma) {
  ifelse(gamma > 1, 1 / gamma, 1)
}

mp_positive_density <- function(x, gamma, sigma2) {
  mp_density(x, gamma, sigma2) / mp_positive_mass(gamma)
}

mp_positive_quantiles <- function(gamma, sigma2, probs, grid_n = 20000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 0) + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_positive_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_fixed_gamma <- function(vals, gamma) {
  vals <- vals[is.finite(vals) & vals > 1e-10]
  empirical_q <- as.numeric(quantile(vals, fit_probs, names = FALSE))

  objective <- function(log_sigma2) {
    sigma2 <- exp(log_sigma2)
    fitted_q <- mp_positive_quantiles(gamma, sigma2, fit_probs, grid_n = grid_n)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- log(unique(pmax(c(median(vals), mean(vals), quantile(vals, 0.75, names = FALSE)), .Machine$double.eps)))
  fits <- lapply(starts, function(start) {
    optim(
      start,
      objective,
      method = "Brent",
      lower = log(min(vals) / 100),
      upper = log(max(vals) * 100)
    )
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

fit_free_gamma <- function(vals, gamma_start) {
  vals <- vals[is.finite(vals) & vals > 1e-10]
  empirical_q <- as.numeric(quantile(vals, fit_probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_positive_quantiles(gamma, sigma2, fit_probs, grid_n = grid_n)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = unique(pmax(c(gamma_start, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8), 1e-4)),
    sigma2 = unique(pmax(c(median(vals), mean(vals), quantile(vals, 0.75, names = FALSE)), .Machine$double.eps))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = log(c(1e-4, min(vals) / 100)),
      upper = log(c(50, max(vals) * 100)),
      control = list(maxit = 1500)
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

fit_bulk_iterative <- function(vals, gamma_theory, fit_type = "free_gamma", max_iter = 4) {
  positive <- vals[is.finite(vals) & vals > 1e-10]
  fit_vals <- positive
  fit <- NULL
  for (iter in seq_len(max_iter)) {
    fit <- if (fit_type == "fixed_gamma") {
      fit_fixed_gamma(fit_vals, gamma_theory)
    } else {
      fit_free_gamma(fit_vals, gamma_theory)
    }
    next_vals <- positive[positive <= fit["mp_upper"]]
    if (length(next_vals) < 50 || length(next_vals) == length(fit_vals)) {
      break
    }
    fit_vals <- next_vals
  }
  c(fit, n_fit = length(fit_vals))
}

eig_df <- read.csv(eig_file)
if ("method" %in% names(eig_df)) {
  raw_df <- eig_df[eig_df$method == "raw_log_normalized", ]
  glm_df <- eig_df[eig_df$method == "glmpca_eta_lowrank", ]
} else {
  raw_df <- eig_df
  raw_df$n_genes <- if ("n_genes" %in% names(raw_df)) raw_df$n_genes else max(raw_df$rank)
  glm_df <- raw_df[0, ]
}
n_values <- sort(unique(raw_df$n_cells))

fit_rows <- list()
curve_rows <- list()
plot_rows <- list()

for (n_cells in n_values) {
  raw_vals <- raw_df$eigenvalue[raw_df$n_cells == n_cells]
  glm_vals <- glm_df$eigenvalue[glm_df$n_cells == n_cells]
  p <- unique(raw_df$n_genes[raw_df$n_cells == n_cells])
  gamma_theory <- p / n_cells

  fixed_fit <- fit_bulk_iterative(raw_vals, gamma_theory, fit_type = "fixed_gamma")
  free_fit <- fit_bulk_iterative(raw_vals, gamma_theory, fit_type = "free_gamma")

  for (fit_type in c("fixed_gamma", "free_gamma")) {
    fit <- if (fit_type == "fixed_gamma") fixed_fit else free_fit
    raw_outliers <- raw_vals[raw_vals > fit["mp_upper"]]
    glm_outliers <- glm_vals[glm_vals > fit["mp_upper"]]

    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      dataset = dataset_name,
      n_cells = n_cells,
      n_genes = p,
      gamma_theory = gamma_theory,
      fit_type = fit_type,
      gamma = fit["gamma"],
      sigma2 = fit["sigma2"],
      mp_lower = fit["mp_lower"],
      mp_upper = fit["mp_upper"],
      objective = fit["objective"],
      convergence = fit["convergence"],
      n_fit = fit["n_fit"],
      n_raw_positive = sum(raw_vals > 1e-10),
      n_raw_spikes = length(raw_outliers),
      raw_spike_excess_sum = sum(raw_outliers - fit["mp_upper"]),
      raw_lambda1_over_mp_upper = max(raw_vals) / fit["mp_upper"],
      n_glm_eta_spikes = length(glm_outliers),
      glm_eta_spike_excess_sum = sum(glm_outliers - fit["mp_upper"]),
      glm_eta_lambda1_over_mp_upper = ifelse(length(glm_vals) > 0, max(glm_vals) / fit["mp_upper"], NA_real_)
    )

    x <- seq(max(fit["mp_lower"], 1e-12), fit["mp_upper"], length.out = 1000)
    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      n_cells = n_cells,
      fit_type = fit_type,
      eigenvalue = x,
      density = mp_positive_density(x, fit["gamma"], fit["sigma2"])
    )
  }

  preferred <- free_fit
  plot_upper <- max(
    preferred["mp_upper"] * 1.25,
    as.numeric(quantile(raw_vals[raw_vals <= preferred["mp_upper"]], 0.995, names = FALSE, na.rm = TRUE))
  )
  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    eigenvalue = raw_vals[raw_vals <= plot_upper]
  )
}

fit_df <- do.call(rbind, fit_rows)
curve_df <- do.call(rbind, curve_rows)
plot_df <- do.call(rbind, plot_rows)

write.csv(fit_df, file.path(out_dir, paste0(dataset_name, "_mp_bulk_fit_spikes_by_n.csv")), row.names = FALSE)
write.csv(curve_df, file.path(out_dir, paste0(dataset_name, "_mp_bulk_fit_curves_by_n.csv")), row.names = FALSE)

free_fit_df <- fit_df[fit_df$fit_type == "free_gamma", ]
free_curve_df <- curve_df[curve_df$fit_type == "free_gamma", ]
free_fit_df$spike_label <- if (nrow(glm_df) > 0) {
  paste0("raw spikes: ", free_fit_df$n_raw_spikes, "\nGLM eta spikes: ", free_fit_df$n_glm_eta_spikes)
} else {
  paste0("raw spikes: ", free_fit_df$n_raw_spikes)
}

plot_df$n_cells_label <- factor(
  plot_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)
free_curve_df$n_cells_label <- factor(
  free_curve_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)
free_fit_df$n_cells_label <- factor(
  free_fit_df$n_cells,
  levels = n_values,
  labels = paste0("n = ", format(n_values, big.mark = ","))
)

theme_set(theme_classic(base_size = 13))

bulk_plot <- ggplot(plot_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_line(data = free_curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.95) +
  geom_vline(data = free_fit_df, aes(xintercept = mp_upper), color = "#C43C39", linetype = "dashed", linewidth = 0.8) +
  geom_text(
    data = free_fit_df,
    aes(
      x = Inf,
      y = Inf,
      label = spike_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.25,
    size = 3.2
  ) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": MP fit to raw log-PCA bulk"),
    subtitle = "Histogram is zoomed to the fitted bulk region; dashed line is fitted MP upper edge",
    x = "Raw log-normalized covariance eigenvalue",
    y = "Density"
  )

spike_counts <- data.frame(n_cells = free_fit_df$n_cells, spike_type = "Raw log-PCA", n_spikes = free_fit_df$n_raw_spikes)
if (nrow(glm_df) > 0) {
  spike_counts <- rbind(
    spike_counts,
    data.frame(n_cells = free_fit_df$n_cells, spike_type = "GLM-PCA eta", n_spikes = free_fit_df$n_glm_eta_spikes)
  )
}

spike_plot <- ggplot(spike_counts, aes(n_cells, n_spikes, color = spike_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = paste0(dataset_name, ": eigenvalues above fitted MP upper edge"),
    x = "Cells",
    y = "Number of spikes",
    color = "Spectrum"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_mp_bulk_fit_by_n.png")), bulk_plot, width = 10, height = 8, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_spike_counts_by_n.png")), spike_plot, width = 8, height = 5, dpi = 220)

print(free_fit_df)
message("Wrote MP bulk fit/spike outputs to ", out_dir)
