#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: compute_snr_proxy_suite.R <eigenvalues_csv> <mp_fit_csv> <out_dir> [dataset_name]", call. = FALSE)
}

eig_file <- args[[1]]
mp_fit_file <- args[[2]]
out_dir <- args[[3]]
dataset_name <- ifelse(length(args) >= 4, args[[4]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

eig_df <- read.csv(eig_file)
if (!("method" %in% names(eig_df))) {
  eig_df$method <- "raw_log_normalized"
}
eig_df <- eig_df[eig_df$method == "raw_log_normalized", ]

mp_df <- read.csv(mp_fit_file)
mp_df <- mp_df[mp_df$fit_type == "free_gamma", ]

effective_rank <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) {
    return(0)
  }
  p <- w / sum(w)
  exp(-sum(p * log(p)))
}

score_one <- function(vals, edge, top_k = 20) {
  vals <- sort(vals[is.finite(vals) & vals > 1e-12], decreasing = TRUE)
  spikes <- vals[vals > edge]
  excess <- pmax(spikes / edge - 1, 0)
  top <- head(vals, min(top_k, length(vals)))

  data.frame(
    n_spikes = length(spikes),
    lambda1_over_edge = vals[1] / edge,
    lambda5_over_edge = vals[min(5, length(vals))] / edge,
    lambda20_over_edge = vals[min(20, length(vals))] / edge,
    top5_mean_over_edge = mean(head(vals, min(5, length(vals))) / edge),
    top20_mean_over_edge = mean(top / edge),
    normalized_spike_excess = sum(excess),
    sqrt_normalized_spike_excess = sqrt(sum(excess)),
    spike_excess_power = sum(spikes - edge),
    spike_excess_power_over_edge = sum(spikes - edge) / edge,
    mean_spike_over_edge = ifelse(length(spikes) > 0, mean(spikes / edge), NA_real_),
    median_spike_over_edge = ifelse(length(spikes) > 0, median(spikes / edge), NA_real_),
    effective_spike_rank = effective_rank(excess),
    total_variance = sum(vals),
    spike_variance_fraction = sum(spikes) / sum(vals),
    spike_excess_fraction = sum(spikes - edge) / sum(vals)
  )
}

rows <- list()
for (i in seq_len(nrow(mp_df))) {
  fit <- mp_df[i, ]
  vals <- eig_df$eigenvalue[eig_df$n_cells == fit$n_cells]
  rows[[length(rows) + 1]] <- cbind(
    dataset = dataset_name,
    n_cells = fit$n_cells,
    n_genes = fit$n_genes,
    mp_gamma = fit$gamma,
    mp_sigma2 = fit$sigma2,
    noise_edge = fit$mp_upper,
    score_one(vals, fit$mp_upper)
  )
}

suite <- do.call(rbind, rows)
rownames(suite) <- NULL
numeric_cols <- setdiff(colnames(suite), "dataset")
for (col in numeric_cols) {
  suite[[col]] <- as.numeric(suite[[col]])
}

out_file <- file.path(out_dir, paste0(dataset_name, "_snr_proxy_suite_by_n.csv"))
write.csv(suite, out_file, row.names = FALSE)

plot_long <- rbind(
  data.frame(n_cells = suite$n_cells, proxy = "n_spikes", value = suite$n_spikes),
  data.frame(n_cells = suite$n_cells, proxy = "lambda1_over_edge", value = suite$lambda1_over_edge),
  data.frame(n_cells = suite$n_cells, proxy = "normalized_spike_excess", value = suite$normalized_spike_excess),
  data.frame(n_cells = suite$n_cells, proxy = "sqrt_normalized_spike_excess", value = suite$sqrt_normalized_spike_excess),
  data.frame(n_cells = suite$n_cells, proxy = "top20_mean_over_edge", value = suite$top20_mean_over_edge),
  data.frame(n_cells = suite$n_cells, proxy = "effective_spike_rank", value = suite$effective_spike_rank),
  data.frame(n_cells = suite$n_cells, proxy = "spike_excess_fraction", value = suite$spike_excess_fraction)
)

plot_long$proxy <- factor(
  plot_long$proxy,
  levels = c(
    "n_spikes",
    "lambda1_over_edge",
    "normalized_spike_excess",
    "sqrt_normalized_spike_excess",
    "top20_mean_over_edge",
    "effective_spike_rank",
    "spike_excess_fraction"
  )
)

theme_set(theme_classic(base_size = 13))

suite_plot <- ggplot(plot_long, aes(n_cells, value)) +
  geom_line(linewidth = 0.85, color = "#1F78B4") +
  geom_point(size = 2, color = "#1F78B4") +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ proxy, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": SNR proxy suite from MP edge"),
    x = "Cells",
    y = "Proxy value"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_snr_proxy_suite_by_n.png")), suite_plot, width = 10, height = 9, dpi = 220)

print(suite)
message("Wrote SNR proxy suite to ", out_file)
