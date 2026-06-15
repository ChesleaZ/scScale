library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/jurkat_all_cells_trimmed_mp_component"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

knee_max_rank <- 100
n_hist_bins <- 120
lower_bins_drop_grid <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 30, 40)
plot_bins_drop <- c(0, 8, 20, 40)
extra_spike_counts <- c(30, 40, 50, 75, 100)

mp_density <- function(x, gamma, sigma2) {
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

mp_quantiles <- function(gamma, sigma2, probs, grid_n = 20000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 1e-12), upper - 1e-12, length.out = grid_n)
  d <- mp_positive_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_free_mp <- function(vals) {
  probs <- seq(0.02, 0.98, by = 0.02)
  empirical_q <- as.numeric(quantile(vals, probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_quantiles(gamma, sigma2, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = c(0.01, 0.03, 0.1, 0.3, 0.7, 1, 1.5, 3, 8, 20),
    sigma2 = unique(as.numeric(c(
      median(vals),
      mean(vals),
      quantile(vals, c(0.25, 0.75), names = FALSE)
    )))
  )

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 2000)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma <- exp(best$par[1])
  sigma2 <- exp(best$par[2])
  list(
    gamma = gamma,
    sigma2 = sigma2,
    lower = sigma2 * (1 - sqrt(gamma))^2,
    upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = best$value
  )
}

find_scree_knee <- function(eigvals, max_rank = 100) {
  eigvals <- sort(eigvals, decreasing = TRUE)
  k <- min(max_rank, length(eigvals) - 1)
  y <- log10(pmax(eigvals[seq_len(k)], .Machine$double.eps))
  x <- seq_along(y)
  x <- (x - min(x)) / (max(x) - min(x))
  y <- (y - min(y)) / (max(y) - min(y))

  x1 <- x[1]
  y1 <- y[1]
  x2 <- tail(x, 1)
  y2 <- tail(y, 1)

  dist <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  which.max(dist)
}

fit_trimmed_component <- function(vals, n_cells, lower_cut, lower_bins_drop, k_spikes) {
  kept <- vals[vals >= lower_cut]
  p_eff <- length(kept)
  count_gamma <- p_eff / n_cells
  fit <- fit_free_mp(kept)

  data.frame(
    lower_bins_drop = lower_bins_drop,
    lower_cut = lower_cut,
    n_for_mp = n_cells,
    p_for_mp = p_eff,
    n_component = p_eff,
    n_lower_trimmed = length(vals) - p_eff,
    gamma_count = count_gamma,
    gamma = fit$gamma,
    sigma2 = fit$sigma2,
    mp_lower = fit$lower,
    mp_upper = fit$upper,
    component_min = min(kept),
    component_median = median(kept),
    component_max = max(kept),
    log_quantile_rmse = sqrt(fit$objective),
    n_above_mp_upper = sum(kept > fit$upper),
    n_below_mp_lower = sum(kept < fit$lower),
    k_spikes = k_spikes
  )
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
summary_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells_summary.csv"))

n_cells <- max(eig_df$n_cells)
all_vals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells], decreasing = TRUE)
p <- length(all_vals)
k_scree <- find_scree_knee(all_vals, max_rank = knee_max_rank)
k_spike_grid <- sort(unique(c(k_scree, extra_spike_counts)))
k_spike_grid <- k_spike_grid[k_spike_grid < p]

fit_list <- list()
component_list <- list()
curve_list <- list()
residual_list <- list()

for (k_spikes in k_spike_grid) {
  residual_vals <- all_vals[(k_spikes + 1):p]
  residual_vals <- residual_vals[residual_vals > 1e-10]
  breaks <- seq(min(residual_vals), max(residual_vals), length.out = n_hist_bins + 1)
  lower_cuts <- breaks[pmin(lower_bins_drop_grid + 1, length(breaks))]

  fit_k <- do.call(rbind, Map(
    function(lower_cut, lower_bins_drop) {
      fit_trimmed_component(residual_vals, n_cells, lower_cut, lower_bins_drop, k_spikes)
    },
    lower_cuts,
    lower_bins_drop_grid
  ))
  fit_k$n_cells <- n_cells
  fit_k$n_genes <- p
  fit_k$n_residual_positive <- length(residual_vals)
  fit_k$n_removed_spikes <- k_spikes
  fit_k$n_removed_nonpositive <- p - k_spikes - length(residual_vals)
  fit_k$n_removed_before_mp <- p - fit_k$p_for_mp
  fit_k$hist_bins <- n_hist_bins
  fit_k$trace_mean_all <- summary_df$mean_eigenvalue[summary_df$n_cells == n_cells]
  fit_list[[as.character(k_spikes)]] <- fit_k

  residual_list[[as.character(k_spikes)]] <- data.frame(k_spikes = k_spikes, eigenvalue = residual_vals)

  component_list[[as.character(k_spikes)]] <- do.call(rbind, lapply(seq_len(nrow(fit_k)), function(i) {
    fit <- fit_k[i, ]
    vals <- residual_vals[residual_vals >= fit$lower_cut]
    data.frame(
      k_spikes = k_spikes,
      lower_bins_drop = fit$lower_bins_drop,
      lower_cut = fit$lower_cut,
      eigenvalue = vals
    )
  }))

  curve_list[[as.character(k_spikes)]] <- do.call(rbind, lapply(seq_len(nrow(fit_k)), function(i) {
    fit <- fit_k[i, ]
    x <- seq(max(fit$mp_lower, 1e-12), fit$mp_upper, length.out = 1000)
    data.frame(
      k_spikes = k_spikes,
      lower_bins_drop = fit$lower_bins_drop,
      lower_cut = fit$lower_cut,
      eigenvalue = x,
      density = mp_positive_density(x, fit$gamma, fit$sigma2)
    )
  }))
}

fit_df <- do.call(rbind, fit_list)
component_df <- do.call(rbind, component_list)
curve_df <- do.call(rbind, curve_list)
residual_df <- do.call(rbind, residual_list)

fit_df$trim_label <- paste0("k=", fit_df$k_spikes, ", drop ", fit_df$lower_bins_drop)
fit_df$k_label <- paste0("k = ", fit_df$k_spikes)
component_df <- merge(component_df, fit_df[, c("k_spikes", "lower_bins_drop", "trim_label", "k_label")], by = c("k_spikes", "lower_bins_drop"))
curve_df <- merge(curve_df, fit_df[, c("k_spikes", "lower_bins_drop", "trim_label", "k_label", "mp_upper")], by = c("k_spikes", "lower_bins_drop"))

best_by_k <- do.call(rbind, lapply(split(fit_df, fit_df$k_spikes), function(df) df[which.min(df$log_quantile_rmse), ]))
best_fit <- fit_df[which.min(fit_df$log_quantile_rmse), ]

plot_fit_df <- subset(fit_df, k_spikes %in% c(k_scree, 50, 100) & lower_bins_drop %in% plot_bins_drop)
plot_component_df <- merge(
  component_df,
  plot_fit_df[, c("k_spikes", "lower_bins_drop")],
  by = c("k_spikes", "lower_bins_drop")
)
plot_curve_df <- merge(
  curve_df,
  plot_fit_df[, c("k_spikes", "lower_bins_drop")],
  by = c("k_spikes", "lower_bins_drop")
)
plot_component_df$trim_label <- factor(plot_component_df$trim_label, levels = plot_fit_df$trim_label)
plot_curve_df$trim_label <- factor(plot_curve_df$trim_label, levels = plot_fit_df$trim_label)
plot_fit_df$trim_label <- factor(plot_fit_df$trim_label, levels = plot_fit_df$trim_label)

write.csv(fit_df, file.path(out_dir, "jurkat_all_cells_trimmed_mp_parameters.csv"), row.names = FALSE)
write.csv(component_df, file.path(out_dir, "jurkat_all_cells_trimmed_components.csv"), row.names = FALSE)

best_fit <- fit_df[which.min(fit_df$log_quantile_rmse), ]
write.csv(best_fit, file.path(out_dir, "jurkat_all_cells_trimmed_mp_best.csv"), row.names = FALSE)
write.csv(best_by_k, file.path(out_dir, "jurkat_all_cells_trimmed_mp_best_by_k.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

sweep_plot <- ggplot(fit_df, aes(lower_bins_drop, log_quantile_rmse, color = k_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  labs(
    title = "All-cell Jurkat trimmed-component MP fit sweep",
    subtitle = paste0("Free empirical MP fit; scree elbow k = ", k_scree, "; also sweeping more aggressive spike removal"),
    x = "Lower residual histogram bins removed",
    y = "Log-quantile RMSE",
    color = "Spikes removed"
  )

component_plot <- ggplot(plot_component_df, aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 90,
    fill = "grey72",
    color = "white",
    linewidth = 0.12
  ) +
  geom_line(
    data = plot_curve_df,
    aes(eigenvalue, density),
    color = "#1F78B4",
    linewidth = 0.9
  ) +
  geom_vline(data = plot_fit_df, aes(xintercept = mp_upper), color = "#C43C39", linetype = "dashed") +
  geom_text(
    data = plot_fit_df,
    aes(
      x = Inf,
      y = Inf,
      label = paste0(
        "n: ", n_component,
        "\nsigma2: ", signif(sigma2, 3),
        "\ngamma: ", signif(gamma, 3),
        "\nrmse: ", signif(log_quantile_rmse, 3)
      )
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.25,
    size = 3.1
  ) +
  facet_wrap(~ trim_label, scales = "free_y", ncol = 2) +
  labs(
    title = "All-cell Jurkat upper residual component vs exact MP",
    subtitle = "Free empirical MP fit; spike count and lower residual-bin trimming swept",
    x = "Covariance eigenvalue",
    y = "Density"
  )

residual_plot <- ggplot(subset(residual_df, k_spikes %in% c(k_scree, 50, 100)), aes(eigenvalue)) +
  geom_histogram(bins = n_hist_bins, fill = "grey75", color = "white", linewidth = 0.12) +
  facet_wrap(~ k_spikes, scales = "free_y", labeller = label_both) +
  labs(
    title = "All-cell Jurkat residual bulk after aggressive spike removal",
    subtitle = paste0("scree elbow k = ", k_scree),
    x = "Covariance eigenvalue after spike removal",
    y = "Count"
  )

ggsave(file.path(out_dir, "jurkat_all_cells_trimmed_mp_sweep.png"), sweep_plot, width = 8, height = 5, dpi = 220)
ggsave(file.path(out_dir, "jurkat_all_cells_trimmed_mp_components.png"), component_plot, width = 11, height = 12, dpi = 220)
ggsave(file.path(out_dir, "jurkat_all_cells_residual_trim_thresholds.png"), residual_plot, width = 8, height = 5, dpi = 220)

print(fit_df)
print(best_fit)
message("Wrote all-cell trimmed MP component outputs to ", out_dir)
