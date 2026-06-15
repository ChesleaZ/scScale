library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"

clip_at <- 2
n_bins <- 220

mp_density_raw <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_positive_quantiles <- function(gamma, sigma2, probs, grid_n = 5000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(lower + 1e-12, upper - 1e-12, length.out = grid_n)
  density <- mp_density_raw(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(density * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_mp <- function(eigvals, clip_at = 2) {
  fit_vals <- eigvals[eigvals > 1e-10 & eigvals <= clip_at]
  probs <- seq(0.02, 0.98, by = 0.02)
  empirical_q <- as.numeric(quantile(fit_vals, probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_positive_quantiles(gamma, sigma2, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = c(0.1, 0.3, 0.7, 1, 1.5, 3, 8),
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
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2

  data.frame(
    mp_fit_gamma = gamma,
    mp_fit_sigma2 = sigma2,
    mp_fit_lower = lower,
    mp_fit_upper = upper,
    mp_fit_objective = best$value,
    n_fit = length(fit_vals)
  )
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
summary_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells_summary.csv"))

cell_counts <- sort(unique(eig_df$n_cells))
labels <- paste0("n = ", format(cell_counts, big.mark = ","))
names(labels) <- cell_counts

fit_df <- do.call(rbind, lapply(cell_counts, function(n_cells) {
  eigvals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  cbind(n_cells = n_cells, fit_mp(eigvals, clip_at = clip_at))
}))

fit_df$mp_fit_density_mass <- ifelse(fit_df$mp_fit_gamma > 1, 1 / fit_df$mp_fit_gamma, 1)
fit_df <- merge(fit_df, summary_df[, c("n_cells", "n_eigenvalues_gt_2")], by = "n_cells")

curve_df <- do.call(rbind, lapply(seq_len(nrow(fit_df)), function(i) {
  row <- fit_df[i, ]
  x <- seq(max(row$mp_fit_lower, 1e-12), min(row$mp_fit_upper, clip_at), length.out = 1000)
  raw_density <- mp_density_raw(x, row$mp_fit_gamma, row$mp_fit_sigma2)
  density <- raw_density / row$mp_fit_density_mass
  data.frame(
    n_cells = row$n_cells,
    eigenvalue = x,
    density = density
  )
}))

eig_df$eigenvalue_clipped <- pmin(eig_df$eigenvalue, clip_at)
eig_df$n_cells_label <- factor(labels[as.character(eig_df$n_cells)], levels = labels)
curve_df$n_cells_label <- factor(labels[as.character(curve_df$n_cells)], levels = labels)
fit_df$n_cells_label <- factor(labels[as.character(fit_df$n_cells)], levels = labels)

write.csv(
  fit_df,
  file.path(out_dir, "jurkat_mp_free_fit_by_cells.csv"),
  row.names = FALSE
)

theme_set(theme_classic(base_size = 13))

facet_plot <- ggplot(eig_df, aes(eigenvalue_clipped)) +
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
  geom_vline(xintercept = clip_at, color = "#C43C39", linetype = "dashed") +
  geom_text(
    data = fit_df,
    aes(x = 1.94, y = Inf, label = paste0(">2: ", n_eigenvalues_gt_2)),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.4,
    size = 3.4
  ) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Jurkat eigenvalue distributions with free MP fits",
    subtitle = paste0(n_bins, " bins; eigenvalues larger than ", clip_at, " are clipped to ", clip_at),
    x = "Covariance eigenvalue, clipped",
    y = "Density"
  )

all_cells <- subset(eig_df, n_cells == max(n_cells))
all_curve <- subset(curve_df, n_cells == max(n_cells))
all_fit <- subset(fit_df, n_cells == max(n_cells))

all_cells_plot <- ggplot(all_cells, aes(eigenvalue_clipped)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = n_bins,
    fill = "grey72",
    color = "white",
    linewidth = 0.12
  ) +
  geom_line(
    data = all_curve,
    aes(eigenvalue, density),
    color = "#1F78B4",
    linewidth = 1
  ) +
  geom_vline(xintercept = clip_at, color = "#C43C39", linetype = "dashed") +
  annotate(
    "text",
    x = 1.94,
    y = Inf,
    label = paste0(">2: ", all_fit$n_eigenvalues_gt_2),
    hjust = 1,
    vjust = 1.4,
    size = 4
  ) +
  labs(
    title = "All-cell Jurkat eigenvalue distribution with free MP fit",
    subtitle = "The empirical noise bulk is visibly not Marchenko-Pastur shaped",
    x = "Covariance eigenvalue, clipped at 2",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_distributions_by_cells_free_mp_fit.png"),
  facet_plot,
  width = 11,
  height = 12,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_all_cells_free_mp_fit.png"),
  all_cells_plot,
  width = 8,
  height = 5,
  dpi = 220
)

print(fit_df)
message("Wrote free MP fit outputs to ", out_dir)
