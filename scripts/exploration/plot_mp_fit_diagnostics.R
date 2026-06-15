#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_mp_fit_diagnostics.R <eigenvalues_csv> <mp_fit_csv> <mp_curves_csv> <out_dir> [dataset_name]", call. = FALSE)
}

eig_file <- args[[1]]
mp_fit_file <- args[[2]]
mp_curves_file <- args[[3]]
out_dir <- args[[4]]
dataset_name <- ifelse(length(args) >= 5, args[[5]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

eig_df <- read.csv(eig_file)
if ("method" %in% names(eig_df)) {
  eig_df <- eig_df[eig_df$method == "raw_log_normalized", ]
}
if (!("n_genes" %in% names(eig_df))) {
  eig_df$n_genes <- max(eig_df$rank)
}

fit_df <- read.csv(mp_fit_file)
fit_df <- fit_df[fit_df$fit_type == "free_gamma", ]
curve_df <- read.csv(mp_curves_file)
curve_df <- curve_df[curve_df$fit_type == "free_gamma", ]

n_values <- sort(unique(eig_df$n_cells))
labels <- paste0("n = ", format(n_values, big.mark = ","))
names(labels) <- n_values

eig_df$n_cells_label <- factor(labels[as.character(eig_df$n_cells)], levels = labels)
fit_df$n_cells_label <- factor(labels[as.character(fit_df$n_cells)], levels = labels)
curve_df$n_cells_label <- factor(labels[as.character(curve_df$n_cells)], levels = labels)

plot_rows <- list()
for (n_cells in n_values) {
  vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  fit <- fit_df[fit_df$n_cells == n_cells, ]
  upper <- fit$mp_upper
  plot_upper <- max(upper * 1.35, quantile(vals[vals <= upper], 0.995, names = FALSE, na.rm = TRUE))
  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    n_cells_label = labels[as.character(n_cells)],
    eigenvalue = vals[vals <= plot_upper]
  )
}
bulk_plot_df <- do.call(rbind, plot_rows)
bulk_plot_df$n_cells_label <- factor(bulk_plot_df$n_cells_label, levels = labels)

theme_set(theme_classic(base_size = 13))

bulk_plot <- ggplot(bulk_plot_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 90, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_line(data = curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.95) +
  geom_vline(data = fit_df, aes(xintercept = mp_upper), color = "#C43C39", linetype = "dashed", linewidth = 0.75) +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": MP bulk diagnostic"),
    subtitle = "Zoomed to bulk region; red dashed line is fitted MP upper edge",
    x = "Covariance eigenvalue",
    y = "Density"
  )

log_plot_df <- eig_df[eig_df$eigenvalue > 1e-10, ]
edge_df <- fit_df
edge_df$log10_mp_upper <- log10(edge_df$mp_upper)

log_hist <- ggplot(log_plot_df, aes(log10(eigenvalue))) +
  geom_histogram(bins = 90, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_vline(data = edge_df, aes(xintercept = log10_mp_upper), color = "#C43C39", linetype = "dashed", linewidth = 0.75) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": full spectrum on log10 eigenvalue scale"),
    subtitle = "Red dashed line is fitted MP upper edge",
    x = "log10 covariance eigenvalue",
    y = "Count"
  )

qq_rows <- list()
cdf_rows <- list()
for (n_cells in n_values) {
  vals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells & eig_df$eigenvalue > 1e-10])
  fit <- fit_df[fit_df$n_cells == n_cells, ]
  curve <- curve_df[curve_df$n_cells == n_cells, ]
  curve <- curve[order(curve$eigenvalue), ]
  dx <- c(diff(curve$eigenvalue), tail(diff(curve$eigenvalue), 1))
  cdf <- cumsum(curve$density * dx)
  cdf <- cdf / max(cdf)

  bulk_vals <- vals[vals <= fit$mp_upper]
  empirical_cdf <- seq_along(bulk_vals) / length(bulk_vals)
  cdf_rows[[length(cdf_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    n_cells_label = labels[as.character(n_cells)],
    eigenvalue = bulk_vals,
    empirical_cdf = empirical_cdf,
    model_cdf = approx(curve$eigenvalue, cdf, xout = bulk_vals, rule = 2)$y
  )

  probs <- seq(0.02, 0.98, length.out = 300)
  qq_rows[[length(qq_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    n_cells_label = labels[as.character(n_cells)],
    prob = probs,
    empirical = as.numeric(quantile(bulk_vals, probs, names = FALSE)),
    model = approx(cdf, curve$eigenvalue, xout = probs, ties = "ordered", rule = 2)$y
  )
}

qq_df <- do.call(rbind, qq_rows)
cdf_df <- do.call(rbind, cdf_rows)
qq_df$n_cells_label <- factor(qq_df$n_cells_label, levels = labels)
cdf_df$n_cells_label <- factor(cdf_df$n_cells_label, levels = labels)

qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
qq_plot <- ggplot(qq_df, aes(model, empirical)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75, color = "#1F78B4") +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": MP bulk Q-Q diagnostic"),
    subtitle = "Only empirical eigenvalues below fitted MP upper edge",
    x = "Fitted MP quantile",
    y = "Empirical bulk quantile"
  )

cdf_plot <- ggplot(cdf_df, aes(eigenvalue)) +
  geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.7) +
  geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.8) +
  facet_wrap(~ n_cells_label, scales = "free_x", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": MP bulk CDF diagnostic"),
    subtitle = "Only empirical eigenvalues below fitted MP upper edge",
    x = "Covariance eigenvalue",
    y = "CDF"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_mp_bulk_zoom_diagnostic.png")), bulk_plot, width = 10, height = 8, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_full_log_spectrum_diagnostic.png")), log_hist, width = 10, height = 8, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_bulk_qq_diagnostic.png")), qq_plot, width = 10, height = 8, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_bulk_cdf_diagnostic.png")), cdf_plot, width = 10, height = 8, dpi = 220)

write.csv(qq_df, file.path(out_dir, paste0(dataset_name, "_mp_bulk_qq_diagnostic.csv")), row.names = FALSE)
write.csv(cdf_df, file.path(out_dir, paste0(dataset_name, "_mp_bulk_cdf_diagnostic.csv")), row.names = FALSE)

message("Wrote MP diagnostics to ", out_dir)
