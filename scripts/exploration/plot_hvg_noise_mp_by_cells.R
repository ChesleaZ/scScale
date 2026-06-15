library(ggplot2)

in_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
out_dir <- "outputs/exploration/jurkat_hvg_noise_mp_by_cells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_bins <- 160
knee_max_rank <- 100

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_quantiles <- function(gamma, sigma2, probs, grid_n = 6000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(lower + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

find_scree_knee <- function(eigvals, max_rank = 100) {
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

fit_mp <- function(vals) {
  fit_vals <- vals[vals > 1e-10]
  probs <- seq(0.02, 0.98, by = 0.02)
  empirical_q <- as.numeric(quantile(fit_vals, probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_quantiles(gamma, sigma2, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = c(0.2, 0.5, 1, 2, 5, 10),
    sigma2 = c(median(fit_vals), mean(fit_vals), quantile(fit_vals, 0.75))
  )

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1500)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma <- exp(best$par[1])
  sigma2 <- exp(best$par[2])

  data.frame(
    mp_gamma = gamma,
    mp_sigma2 = sigma2,
    mp_lower = sigma2 * (1 - sqrt(gamma))^2,
    mp_upper = sigma2 * (1 + sqrt(gamma))^2,
    mp_objective = best$value,
    n_positive_fit = length(fit_vals)
  )
}

eig_df <- read.csv(in_file)
cell_counts <- sort(unique(eig_df$n_cells))

noise_list <- vector("list", length(cell_counts))
summary_list <- vector("list", length(cell_counts))
curve_list <- vector("list", length(cell_counts))

for (i in seq_along(cell_counts)) {
  n_cells <- cell_counts[i]
  eigvals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  eigvals <- sort(eigvals, decreasing = TRUE)

  spike_count <- find_scree_knee(eigvals, max_rank = knee_max_rank)
  noise_upper <- eigvals[spike_count + 1]
  noise_vals <- eigvals[(spike_count + 1):length(eigvals)]

  fit <- fit_mp(noise_vals)

  x_max <- min(max(noise_vals), max(fit$mp_upper, noise_upper))
  x <- seq(max(fit$mp_lower, 1e-12), x_max, length.out = 1000)
  positive_mass <- mean(noise_vals > 1e-10)

  curve_list[[i]] <- data.frame(
    n_cells = n_cells,
    eigenvalue = x,
    density = mp_density(x, fit$mp_gamma, fit$mp_sigma2) * positive_mass
  )

  noise_list[[i]] <- data.frame(
    n_cells = n_cells,
    rank = seq_along(noise_vals) + spike_count,
    eigenvalue = noise_vals
  )

  summary_list[[i]] <- cbind(
    data.frame(
      n_cells = n_cells,
      spike_count = spike_count,
      noise_upper = noise_upper,
      n_noise = length(noise_vals),
      n_noise_zero = sum(noise_vals <= 1e-10)
    ),
    fit
  )
}

noise_df <- do.call(rbind, noise_list)
summary_df <- do.call(rbind, summary_list)
curve_df <- do.call(rbind, curve_list)

labels <- paste0("n = ", format(cell_counts, big.mark = ","))
names(labels) <- cell_counts

noise_df$n_cells_label <- factor(labels[as.character(noise_df$n_cells)], levels = labels)
curve_df$n_cells_label <- factor(labels[as.character(curve_df$n_cells)], levels = labels)
summary_df$n_cells_label <- factor(labels[as.character(summary_df$n_cells)], levels = labels)

write.csv(
  noise_df,
  file.path(out_dir, "jurkat_hvg_noise_eigenvalues_by_cells.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(out_dir, "jurkat_hvg_noise_threshold_mp_fit_by_cells.csv"),
  row.names = FALSE
)

theme_set(theme_classic(base_size = 13))

hist_plot <- ggplot(noise_df, aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = n_bins,
    fill = "grey72",
    color = "white",
    linewidth = 0.12
  ) +
  geom_line(
    data = curve_df,
    aes(eigenvalue, density),
    color = "#1F78B4",
    linewidth = 0.9
  ) +
  geom_vline(
    data = summary_df,
    aes(xintercept = noise_upper),
    color = "#C43C39",
    linetype = "dashed"
  ) +
  geom_text(
    data = summary_df,
    aes(
      x = Inf,
      y = Inf,
      label = paste0("spikes: ", spike_count, "\nthreshold: ", signif(noise_upper, 3))
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.25,
    size = 3.2
  ) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  coord_cartesian(xlim = c(0, max(summary_df$noise_upper) * 1.05)) +
  labs(
    title = "HVG Jurkat noise eigenvalue distributions with free MP overlays",
    subtitle = paste0("Noise threshold is selected automatically from the scree elbow among the top ", knee_max_rank, " ranks"),
    x = "Covariance eigenvalue",
    y = "Density"
  )

threshold_plot <- ggplot(summary_df, aes(n_cells, spike_count)) +
  geom_line(linewidth = 0.9, color = "grey25") +
  geom_point(size = 2.2, color = "#1F78B4") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Automatically selected spike count",
    x = "Number of cells",
    y = "Spike count"
  )

ggsave(
  file.path(out_dir, "jurkat_hvg_noise_distributions_by_cells_free_mp.png"),
  hist_plot,
  width = 11,
  height = 12,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_hvg_noise_spike_count_by_cells.png"),
  threshold_plot,
  width = 8,
  height = 5,
  dpi = 220
)

print(summary_df)
message("Wrote outputs to ", out_dir)
