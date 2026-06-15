#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: compute_mp_snr_proxy_from_fit.R <eigenvalues_csv> <mp_fit_csv> <out_dir> [dataset_name]", call. = FALSE)
}

eig_file <- args[[1]]
mp_fit_file <- args[[2]]
out_dir <- args[[3]]
dataset_name <- ifelse(length(args) >= 4, args[[4]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

eig_df <- read.csv(eig_file)
mp_df <- read.csv(mp_fit_file)
mp_df <- mp_df[mp_df$fit_type == "free_gamma", ]
if (!("method" %in% names(eig_df))) {
  eig_df$method <- "raw_log_normalized"
}

score_spectrum <- function(vals, edge) {
  vals <- vals[is.finite(vals) & vals > 1e-12]
  spikes <- vals[vals > edge]
  data.frame(
    n_spikes = length(spikes),
    lambda1 = max(vals),
    lambda1_over_edge = max(vals) / edge,
    spike_excess_sum = sum(spikes - edge),
    normalized_spike_excess = sum(spikes / edge - 1),
    mean_spike_over_edge = ifelse(length(spikes) > 0, mean(spikes / edge), NA_real_),
    median_spike_over_edge = ifelse(length(spikes) > 0, median(spikes / edge), NA_real_)
  )
}

rows <- list()
for (i in seq_len(nrow(mp_df))) {
  fit <- mp_df[i, ]
  n_cells <- fit$n_cells
  edge <- fit$mp_upper

  for (method in intersect(c("raw_log_normalized", "glmpca_eta_lowrank"), unique(eig_df$method))) {
    vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells & eig_df$method == method]
    score <- score_spectrum(vals, edge)
    rows[[length(rows) + 1]] <- cbind(
      dataset = dataset_name,
      n_cells = n_cells,
      n_genes = fit$n_genes,
      spectrum = method,
      mp_gamma = fit$gamma,
      mp_sigma2 = fit$sigma2,
      noise_eigenvalue_mp_upper = edge,
      mp_lower = fit$mp_lower,
      score
    )
  }
}

proxy_df <- do.call(rbind, rows)
rownames(proxy_df) <- NULL
numeric_cols <- setdiff(colnames(proxy_df), c("dataset", "spectrum"))
for (col in numeric_cols) {
  proxy_df[[col]] <- as.numeric(proxy_df[[col]])
}

out_file <- file.path(out_dir, paste0(dataset_name, "_mp_snr_proxy_by_n.csv"))
write.csv(proxy_df, out_file, row.names = FALSE)

plot_df <- proxy_df
plot_df$spectrum_label <- factor(
  plot_df$spectrum,
  levels = c("raw_log_normalized", "glmpca_eta_lowrank"),
  labels = c("Raw log-PCA", "GLM-PCA eta")
)

theme_set(theme_classic(base_size = 13))

proxy_plot <- ggplot(plot_df, aes(n_cells, normalized_spike_excess, color = spectrum_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": MP-edge normalized spike excess"),
    subtitle = "sum(max(lambda / MP upper edge - 1, 0)); MP edge is the noise eigenvalue proxy",
    x = "Cells",
    y = "Normalized spike excess (log scale)",
    color = "Spectrum"
  )

lambda_plot <- ggplot(plot_df, aes(n_cells, lambda1_over_edge, color = spectrum_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": leading spike over MP edge"),
    x = "Cells",
    y = "lambda1 / MP upper edge (log scale)",
    color = "Spectrum"
  )

edge_plot <- ggplot(mp_df, aes(n_cells, mp_upper)) +
  geom_line(linewidth = 0.9, color = "#1F78B4") +
  geom_point(size = 2.2, color = "#1F78B4") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": fitted MP upper edge"),
    subtitle = "Noise eigenvalue proxy from raw log-PCA bulk",
    x = "Cells",
    y = "MP upper edge (log scale)"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_mp_normalized_spike_excess_by_n.png")), proxy_plot, width = 8.5, height = 5.25, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_lambda1_over_edge_by_n.png")), lambda_plot, width = 8.5, height = 5.25, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_mp_noise_edge_by_n.png")), edge_plot, width = 8, height = 5, dpi = 220)

print(proxy_df)
message("Wrote MP SNR proxy table to ", out_file)
