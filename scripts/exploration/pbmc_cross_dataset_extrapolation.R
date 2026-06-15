#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- ifelse(length(args) >= 1, args[[1]], "outputs/exploration/pbmc_cross_dataset_extrapolation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

pbmc_inputs <- data.frame(
  dataset_id = c(
    "10x_PBMC_Multiome_3k",
    "SmartSeq3_PBMC_UMI",
    "FedscGen_PBMC",
    "Parse_1M_PBMC_T1D"
  ),
  dataset_label = c(
    "10x PBMC Multiome 3k",
    "Smart-seq3 PBMC",
    "FedscGen PBMC",
    "Parse 1M PBMC T1D"
  ),
  plot_label = c("10x 3k", "SS3", "FedscGen", "Parse"),
  report_dir = file.path(
    "reports/datasets",
    c(
      "10x_PBMC_Multiome_3k",
      "SmartSeq3_PBMC_UMI",
      "FedscGen_PBMC",
      "Parse_1M_PBMC_T1D"
    )
  )
)

effective_rank <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) {
    return(0)
  }
  p <- w / sum(w)
  exp(-sum(p * log(p)))
}

score_dataset <- function(dataset_id, dataset_label, plot_label, report_dir) {
  eig_path <- file.path(report_dir, "pca_eigenvalues.csv")
  mp_path <- file.path(report_dir, "mp_fit_parameters.csv")
  if (!file.exists(eig_path) || !file.exists(mp_path)) {
    warning("Skipping ", dataset_id, ": missing PCA eigenvalues or MP fit parameters")
    return(NULL)
  }

  eig <- read.csv(eig_path)
  mp <- read.csv(mp_path)
  vals <- sort(eig$eigenvalue[is.finite(eig$eigenvalue) & eig$eigenvalue > 1e-12], decreasing = TRUE)
  edge <- mp$mp_upper[1]
  spikes <- vals[vals > edge]
  excess <- pmax(spikes / edge - 1, 0)
  top20 <- head(vals, min(20, length(vals)))

  data.frame(
    dataset_id = dataset_id,
    dataset_label = dataset_label,
    plot_label = plot_label,
    n_cells = mp$n_cells[1],
    n_eigenvalues = mp$n_eigenvalues[1],
    mp_gamma = mp$gamma[1],
    mp_sigma2 = mp$sigma2[1],
    noise_edge = edge,
    n_spikes = length(spikes),
    report_n_outliers_above_mp = mp$n_outliers_above_mp[1],
    lambda1_over_edge = vals[1] / edge,
    top20_mean_over_edge = mean(top20 / edge),
    normalized_spike_excess = sum(excess),
    sqrt_normalized_spike_excess = sqrt(sum(excess)),
    spike_excess_power = sum(spikes - edge),
    spike_excess_power_over_edge = sum(spikes - edge) / edge,
    spike_variance_fraction = sum(spikes) / sum(vals),
    spike_excess_fraction = sum(spikes - edge) / sum(vals),
    effective_spike_rank = effective_rank(excess),
    total_variance = sum(vals)
  )
}

metric_specs <- data.frame(
  metric = c(
    "lambda1_over_edge",
    "top20_mean_over_edge",
    "sqrt_normalized_spike_excess",
    "n_spikes",
    "noise_edge",
    "spike_excess_fraction",
    "effective_spike_rank",
    "spike_variance_fraction"
  ),
  interpretation = c(
    "Top eigenvalue relative to MP edge",
    "Mean top-20 eigenvalue relative to MP edge",
    "Square-root total spike excess above MP edge",
    "Number of MP-edge spikes",
    "Fitted MP upper edge",
    "Spike excess as fraction of total variance",
    "Effective rank of spike excess",
    "Spike variance fraction"
  )
)

pbmc_metrics <- do.call(rbind, lapply(seq_len(nrow(pbmc_inputs)), function(i) {
  score_dataset(pbmc_inputs$dataset_id[i], pbmc_inputs$dataset_label[i], pbmc_inputs$plot_label[i], pbmc_inputs$report_dir[i])
}))
pbmc_metrics <- pbmc_metrics[order(pbmc_metrics$n_cells), ]

fit_power_law <- function(train_df, target_df, metric) {
  train <- train_df[is.finite(train_df[[metric]]) & train_df[[metric]] > 0, ]
  if (nrow(train) < 2 || !is.finite(target_df[[metric]]) || target_df[[metric]] <= 0) {
    return(NULL)
  }

  fit_df <- data.frame(log_value = log10(train[[metric]]), log_n = log10(train$n_cells))
  fit <- lm(log_value ~ log_n, data = fit_df)
  new_df <- data.frame(log_n = log10(target_df$n_cells))
  if (df.residual(fit) > 0) {
    pred <- predict(fit, newdata = new_df, interval = "prediction", level = 0.80)
    predicted <- 10^pred[, "fit"]
    low80 <- 10^pred[, "lwr"]
    high80 <- 10^pred[, "upr"]
  } else {
    predicted <- 10^predict(fit, newdata = new_df)
    low80 <- NA_real_
    high80 <- NA_real_
  }

  data.frame(
    target_dataset_id = target_df$dataset_id,
    target_dataset_label = target_df$dataset_label,
    target_n_cells = target_df$n_cells,
    train_dataset_ids = paste(train$dataset_id, collapse = ";"),
    train_dataset_labels = paste(train$dataset_label, collapse = "; "),
    n_train_datasets = nrow(train),
    train_n_min = min(train$n_cells),
    train_n_max = max(train$n_cells),
    metric = metric,
    observed = target_df[[metric]],
    predicted = predicted,
    pred_low80 = low80,
    pred_high80 = high80,
    abs_pct_error = abs(predicted / target_df[[metric]] - 1),
    log10_error = log10(predicted) - log10(target_df[[metric]]),
    exponent = unname(coef(fit)[["log_n"]]),
    train_r2 = summary(fit)$r.squared,
    model = "smaller_pbmc_power_law"
  )
}

transfer_rows <- list()
for (target_i in seq_len(nrow(pbmc_metrics))) {
  target <- pbmc_metrics[target_i, ]
  smaller <- pbmc_metrics[pbmc_metrics$n_cells < target$n_cells, ]
  for (metric in metric_specs$metric) {
    transfer_rows[[length(transfer_rows) + 1]] <- fit_power_law(smaller, target, metric)
  }
}
transfer <- do.call(rbind, transfer_rows[!vapply(transfer_rows, is.null, logical(1))])
transfer <- merge(transfer, metric_specs, by = "metric", all.x = TRUE)
transfer <- transfer[order(transfer$target_n_cells, transfer$metric), ]

summary_by_target <- do.call(rbind, lapply(split(transfer, transfer$target_dataset_label), function(d) {
  data.frame(
    target_dataset_label = d$target_dataset_label[1],
    target_n_cells = d$target_n_cells[1],
    n_metrics = nrow(d),
    mean_abs_pct_error = mean(d$abs_pct_error),
    median_abs_pct_error = median(d$abs_pct_error),
    max_abs_pct_error = max(d$abs_pct_error),
    mean_abs_log10_error = mean(abs(d$log10_error))
  )
}))
summary_by_target <- summary_by_target[order(summary_by_target$target_n_cells), ]

metric_summary <- do.call(rbind, lapply(split(transfer, transfer$metric), function(d) {
  data.frame(
    metric = d$metric[1],
    interpretation = d$interpretation[1],
    n_targets = nrow(d),
    mean_abs_pct_error = mean(d$abs_pct_error),
    median_abs_pct_error = median(d$abs_pct_error),
    max_abs_pct_error = max(d$abs_pct_error),
    mean_abs_log10_error = mean(abs(d$log10_error))
  )
}))
metric_summary <- metric_summary[order(metric_summary$mean_abs_pct_error), ]

next_n_cells <- as.numeric(strsplit(Sys.getenv("PBMC_CROSS_DATASET_PREDICT_N_CELLS", "500000,1000000"), ",")[[1]])
forecast_rows <- list()
for (metric in metric_specs$metric) {
  train <- pbmc_metrics[is.finite(pbmc_metrics[[metric]]) & pbmc_metrics[[metric]] > 0, ]
  if (nrow(train) < 2) {
    next
  }
  fit_df <- data.frame(log_value = log10(train[[metric]]), log_n = log10(train$n_cells))
  fit <- lm(log_value ~ log_n, data = fit_df)
  future <- data.frame(n_cells = next_n_cells[next_n_cells > max(train$n_cells)])
  if (nrow(future) == 0) {
    next
  }
  new_df <- data.frame(log_n = log10(future$n_cells))
  if (df.residual(fit) > 0) {
    pred <- predict(fit, newdata = new_df, interval = "prediction", level = 0.80)
    predicted <- 10^pred[, "fit"]
    low80 <- 10^pred[, "lwr"]
    high80 <- 10^pred[, "upr"]
  } else {
    predicted <- 10^predict(fit, newdata = new_df)
    low80 <- NA_real_
    high80 <- NA_real_
  }
  forecast_rows[[length(forecast_rows) + 1]] <- data.frame(
    metric = metric,
    fit_dataset_ids = paste(train$dataset_id, collapse = ";"),
    fit_n_min = min(train$n_cells),
    fit_n_max = max(train$n_cells),
    n_fit_datasets = nrow(train),
    target_n_cells = future$n_cells,
    predicted = predicted,
    pred_low80 = low80,
    pred_high80 = high80,
    exponent = unname(coef(fit)[["log_n"]]),
    train_r2 = summary(fit)$r.squared
  )
}
forecast <- do.call(rbind, forecast_rows)
forecast <- merge(forecast, metric_specs, by = "metric", all.x = TRUE)

write.csv(pbmc_metrics, file.path(out_dir, "pbmc_natural_dataset_metrics.csv"), row.names = FALSE)
write.csv(transfer, file.path(out_dir, "pbmc_cross_dataset_transfer_validation.csv"), row.names = FALSE)
write.csv(summary_by_target, file.path(out_dir, "pbmc_cross_dataset_target_summary.csv"), row.names = FALSE)
write.csv(metric_summary, file.path(out_dir, "pbmc_cross_dataset_metric_summary.csv"), row.names = FALSE)
write.csv(forecast, file.path(out_dir, "pbmc_cross_dataset_future_predictions.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 12))

long_metrics <- do.call(rbind, lapply(metric_specs$metric, function(metric) {
  data.frame(
    dataset_label = pbmc_metrics$dataset_label,
    plot_label = pbmc_metrics$plot_label,
    n_cells = pbmc_metrics$n_cells,
    metric = metric,
    value = pbmc_metrics[[metric]],
    interpretation = metric_specs$interpretation[metric_specs$metric == metric]
  )
}))

short_metric_labels <- c(
  lambda1_over_edge = "lambda1 / edge",
  top20_mean_over_edge = "top20 mean / edge",
  sqrt_normalized_spike_excess = "sqrt spike excess",
  n_spikes = "MP-edge spikes",
  noise_edge = "MP upper edge",
  spike_excess_fraction = "spike excess frac",
  effective_spike_rank = "effective spike rank",
  spike_variance_fraction = "spike variance frac"
)
long_metrics$metric_label <- factor(
  unname(short_metric_labels[long_metrics$metric]),
  levels = unname(short_metric_labels[metric_specs$metric])
)

metrics_plot <- ggplot(long_metrics, aes(n_cells, value)) +
  geom_line(color = "#1F78B4", linewidth = 0.7) +
  geom_point(size = 2.3, color = "#222222") +
  geom_text(aes(label = plot_label), hjust = -0.08, vjust = 0.5, size = 2.8, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma, limits = range(long_metrics$n_cells) * c(0.86, 1.45)) +
  scale_y_log10() +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Natural PBMC datasets: spectrum metrics by cell count",
    subtitle = "Each point is a different PBMC dataset/report, not a subsample.",
    x = "Cells",
    y = "Metric value"
  )

ggsave(file.path(out_dir, "pbmc_natural_dataset_scaling.png"), metrics_plot, width = 11, height = 9, dpi = 220)

transfer_plot <- ggplot(transfer, aes(reorder(interpretation, abs_pct_error), abs_pct_error, fill = target_dataset_label)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_y_log10(labels = scales::percent) +
  labs(
    title = "PBMC-to-PBMC transfer error without within-dataset subsampling",
    subtitle = "For each target, fit uses only smaller natural PBMC datasets. X-axis is log-scaled because FedscGen has a large two-source extrapolation failure.",
    x = NULL,
    y = "Absolute percent error",
    fill = "Target"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "pbmc_cross_dataset_transfer_error.png"), transfer_plot, width = 11, height = 7, dpi = 220)

print(pbmc_metrics)
print(summary_by_target)
print(metric_summary)
message("Wrote PBMC cross-dataset extrapolation outputs to ", out_dir)
