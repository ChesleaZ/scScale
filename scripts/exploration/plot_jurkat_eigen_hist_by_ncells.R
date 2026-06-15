library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  eig_file <- args[[1]]
  out_dir <- args[[2]]
  dataset_label <- args[[3]]
} else {
  eig_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
  out_dir <- "outputs/exploration/jurkat_mp_by_ncells"
  dataset_label <- "Jurkat"
}
mp_file <- file.path(out_dir, "jurkat_mp_by_ncells_curves.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eig_df <- read.csv(eig_file)
mp_curves <- read.csv(mp_file)
mp_params <- read.csv(file.path(out_dir, "jurkat_mp_by_ncells_parameters.csv"))

eig_df <- eig_df[is.finite(eig_df$eigenvalue), ]
eig_df$positive <- eig_df$eigenvalue > 1e-10
n_values <- sort(unique(eig_df$n_cells))
labels <- paste0("n = ", format(n_values, big.mark = ","))

eig_df$n_cells_label <- factor(eig_df$n_cells, levels = n_values, labels = labels)
mp_curves$n_cells_label <- factor(mp_curves$n_cells, levels = n_values, labels = labels)

positive_df <- eig_df[eig_df$positive, ]
free_curves <- mp_curves[mp_curves$fit_type == "free_gamma", ]

theme_set(theme_classic(base_size = 13))

raw_hist <- ggplot(positive_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 140, fill = "grey75", color = "white", linewidth = 0.1) +
  geom_line(data = free_curves, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.85) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_label, " positive eigenvalue histograms across cell counts"),
    subtitle = "Blue curve is free-gamma MP fit",
    x = "Covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "jurkat_positive_eigen_hist_free_mp_by_ncells.png"), raw_hist, width = 10, height = 8, dpi = 220)

zoom_limits <- aggregate(eigenvalue ~ n_cells + n_cells_label, positive_df, function(x) {
  as.numeric(quantile(x, 0.98, names = FALSE))
})
names(zoom_limits)[3] <- "xmax"
positive_zoom <- merge(positive_df, zoom_limits, by = c("n_cells", "n_cells_label"))
positive_zoom <- positive_zoom[positive_zoom$eigenvalue <= positive_zoom$xmax, ]
free_zoom <- merge(free_curves, zoom_limits, by = c("n_cells", "n_cells_label"))
free_zoom <- free_zoom[free_zoom$eigenvalue <= free_zoom$xmax, ]

zoom_hist <- ggplot(positive_zoom, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey75", color = "white", linewidth = 0.1) +
  geom_line(data = free_zoom, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.85) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_label, " positive eigenvalue histograms, lower 98%"),
    subtitle = "Zoom removes the largest positive eigenvalues so the small-eigenvalue pileup is visible",
    x = "Covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "jurkat_positive_eigen_hist_lower98_free_mp_by_ncells.png"), zoom_hist, width = 10, height = 8, dpi = 220)

log_hist <- ggplot(positive_df, aes(log10(eigenvalue))) +
  geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey75", color = "white", linewidth = 0.1) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_label, " log10 positive eigenvalue histograms across cell counts"),
    x = "log10 covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "jurkat_log10_positive_eigen_hist_by_ncells.png"), log_hist, width = 10, height = 8, dpi = 220)

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

free_params <- mp_params[mp_params$fit_type == "free_gamma", ]
log_curve <- do.call(rbind, lapply(n_values, function(n_cells) {
  vals <- positive_df$eigenvalue[positive_df$n_cells == n_cells]
  y <- seq(
    as.numeric(quantile(log10(vals), 0.001, names = FALSE)),
    as.numeric(quantile(log10(vals), 0.999, names = FALSE)),
    length.out = 4000
  )
  x <- 10^y
  fit <- free_params[free_params$n_cells == n_cells, ]
  data.frame(
    n_cells = n_cells,
    n_cells_label = factor(
      n_cells,
      levels = n_values,
      labels = labels
    ),
    log10_eigenvalue = y,
    density = mp_density(x, fit$gamma, fit$sigma2) * x * log(10)
  )
}))

log_hist_free_mp <- ggplot(positive_df, aes(log10(eigenvalue))) +
  geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey75", color = "white", linewidth = 0.1) +
  geom_line(data = log_curve, aes(log10_eigenvalue, density), color = "#1F78B4", linewidth = 0.85) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_label, " log10 eigenvalue histograms with transformed free-gamma MP density"),
    subtitle = "Blue density includes Jacobian: f_log10(x) = f_MP(10^x) * 10^x * log(10)",
    x = "log10 covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "jurkat_log10_positive_eigen_hist_free_mp_overlay_by_ncells.png"), log_hist_free_mp, width = 10, height = 8, dpi = 220)

write.csv(log_curve, file.path(out_dir, "jurkat_free_gamma_mp_log10_density_curves.csv"), row.names = FALSE)

near_zero_summary <- do.call(rbind, lapply(n_values, function(n_cells) {
  vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  pos <- vals[vals > 1e-10]
  data.frame(
    n_cells = n_cells,
    p = length(vals),
    n_zero = sum(vals <= 1e-10),
    zero_fraction = mean(vals <= 1e-10),
    q001 = as.numeric(quantile(pos, 0.001, names = FALSE)),
    q005 = as.numeric(quantile(pos, 0.005, names = FALSE)),
    q01 = as.numeric(quantile(pos, 0.01, names = FALSE)),
    q05 = as.numeric(quantile(pos, 0.05, names = FALSE)),
    q10 = as.numeric(quantile(pos, 0.10, names = FALSE)),
    median_positive = median(pos),
    mean_positive = mean(pos)
  )
}))

write.csv(near_zero_summary, file.path(out_dir, "jurkat_near_zero_eigen_summary_by_ncells.csv"), row.names = FALSE)

print(near_zero_summary)
message("Wrote histogram diagnostics to ", out_dir)
