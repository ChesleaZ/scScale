#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- ifelse(length(args) >= 1, args[[1]], "outputs/exploration/consortium_scaling_extrapolation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
next_n_cells <- as.numeric(strsplit(Sys.getenv("SCALING_PREDICT_N_CELLS", "500000,1000000"), ",")[[1]])

library(ggplot2)

read_existing <- function(path, label, source_type) {
  if (!file.exists(path)) {
    warning("Missing input: ", path)
    return(NULL)
  }
  df <- read.csv(path)
  df$analysis_dataset <- label
  df$source_type <- source_type
  df
}

inputs <- list(
  read_existing(
    "outputs/exploration/jurkat_log_pca_snr_by_n/jurkat_log_pca_snr_proxy_suite_by_n.csv",
    "Jurkat log-PCA",
    "small_to_full_validation"
  ),
  read_existing(
    "outputs/exploration/parse_pbmc_log_pca_200k_snr/parse_pbmc_log_pca_200k_snr_proxy_suite_by_n.csv",
    "Parse PBMC log-PCA",
    "small_to_large_validation"
  ),
  read_existing(
    "outputs/exploration/parse_pbmc_glmpca_snr_by_n/parse_pbmc_glmpca_spike_bulk_proxy_by_n.csv",
    "Parse PBMC GLM-PCA",
    "small_n_probe"
  )
)
inputs <- inputs[!vapply(inputs, is.null, logical(1))]
all_cols <- unique(unlist(lapply(inputs, names)))
inputs <- lapply(inputs, function(df) {
  missing_cols <- setdiff(all_cols, names(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  df[, all_cols]
})
raw <- do.call(rbind, inputs)

metric_specs <- data.frame(
  metric = c(
    "lambda1_over_edge",
    "top20_mean_over_edge",
    "sqrt_normalized_spike_excess",
    "n_spikes",
    "noise_edge",
    "spike_excess_fraction",
    "raw_log_lambda1_over_median",
    "glm_eta_lambda1_over_raw_median",
    "glm_eta_lambda1_over_eta_nonzero_median"
  ),
  interpretation = c(
    "Top eigenvalue relative to MP edge",
    "Mean top-20 eigenvalue relative to MP edge",
    "Square-root total spike excess above MP edge",
    "Number of MP-edge spikes",
    "Fitted MP upper edge",
    "Spike excess as fraction of total variance",
    "Raw log-PCA top eigenvalue over bulk median",
    "GLM-PCA eta top eigenvalue over raw bulk median",
    "GLM-PCA eta top eigenvalue over eta nonzero median"
  )
)

fit_prefix <- function(df, metric, train_n_max = Inf, min_train = 3) {
  keep <- is.finite(df[[metric]]) & df[[metric]] > 0 & is.finite(df$n_cells) & df$n_cells > 0
  d <- df[keep, c("analysis_dataset", "source_type", "n_cells", metric)]
  names(d)[names(d) == metric] <- "observed"
  d <- d[order(d$n_cells), ]
  if (nrow(d) < min_train + 1) {
    return(NULL)
  }

  train <- d[d$n_cells <= train_n_max, ]
  if (nrow(train) < min_train) {
    train <- d[seq_len(min_train), ]
  }
  if (nrow(train) >= nrow(d)) {
    train <- d[seq_len(max(min_train, floor(nrow(d) / 2))), ]
  }
  test <- d[d$n_cells > max(train$n_cells), ]
  if (nrow(test) == 0) {
    return(NULL)
  }

  model <- lm(log10(observed) ~ log10(n_cells), data = train)
  pred_all <- predict(model, newdata = d, interval = "prediction", level = 0.80)
  d$predicted <- 10^pred_all[, "fit"]
  d$pred_low80 <- 10^pred_all[, "lwr"]
  d$pred_high80 <- 10^pred_all[, "upr"]
  d$split <- ifelse(d$n_cells <= max(train$n_cells), "fit", "heldout")
  d$metric <- metric
  d$train_n_max <- max(train$n_cells)
  d$exponent <- unname(coef(model)[["log10(n_cells)"]])
  d$intercept <- unname(coef(model)[["(Intercept)"]])
  d$train_r2 <- summary(model)$r.squared
  d$ape <- abs(d$predicted / d$observed - 1)
  d$log10_error <- log10(d$predicted) - log10(d$observed)
  d
}

choose_train_max <- function(n_cells) {
  n <- sort(unique(n_cells[is.finite(n_cells)]))
  if (length(n) < 4) {
    return(Inf)
  }
  n[max(3, floor(length(n) / 2))]
}

predictions <- list()
for (dataset_name in unique(raw$analysis_dataset)) {
  df <- raw[raw$analysis_dataset == dataset_name, ]
  train_n_max <- choose_train_max(df$n_cells)
  for (metric in metric_specs$metric) {
    if (metric %in% names(df)) {
      predictions[[length(predictions) + 1]] <- fit_prefix(df, metric, train_n_max)
    }
  }
}
pred <- do.call(rbind, predictions[!vapply(predictions, is.null, logical(1))])
pred <- merge(pred, metric_specs, by = "metric", all.x = TRUE)

fit_all_and_predict <- function(df, metric, future_n) {
  keep <- is.finite(df[[metric]]) & df[[metric]] > 0 & is.finite(df$n_cells) & df$n_cells > 0
  d <- df[keep, c("analysis_dataset", "source_type", "n_cells", metric)]
  names(d)[names(d) == metric] <- "observed"
  d <- d[order(d$n_cells), ]
  if (nrow(d) < 3) {
    return(NULL)
  }

  model <- lm(log10(observed) ~ log10(n_cells), data = d)
  future <- data.frame(n_cells = future_n[future_n > max(d$n_cells)])
  if (nrow(future) == 0) {
    return(NULL)
  }
  pred_int <- predict(model, newdata = future, interval = "prediction", level = 0.80)
  data.frame(
    analysis_dataset = d$analysis_dataset[1],
    source_type = d$source_type[1],
    metric = metric,
    fit_n_min = min(d$n_cells),
    fit_n_max = max(d$n_cells),
    n_fit_points = nrow(d),
    target_n_cells = future$n_cells,
    predicted = 10^pred_int[, "fit"],
    pred_low80 = 10^pred_int[, "lwr"],
    pred_high80 = 10^pred_int[, "upr"],
    exponent = unname(coef(model)[["log10(n_cells)"]]),
    train_r2 = summary(model)$r.squared
  )
}

future_predictions <- list()
for (dataset_name in unique(raw$analysis_dataset)) {
  df <- raw[raw$analysis_dataset == dataset_name, ]
  for (metric in metric_specs$metric) {
    if (metric %in% names(df)) {
      future_predictions[[length(future_predictions) + 1]] <- fit_all_and_predict(df, metric, next_n_cells)
    }
  }
}
future_pred <- do.call(rbind, future_predictions[!vapply(future_predictions, is.null, logical(1))])
future_pred <- merge(future_pred, metric_specs, by = "metric", all.x = TRUE)

heldout <- pred[pred$split == "heldout", ]
summary_one <- function(d) {
  data.frame(
    mean_abs_pct_error = mean(d$ape),
    median_abs_pct_error = median(d$ape),
    max_abs_pct_error = max(d$ape),
    mean_abs_log10_error = mean(abs(d$log10_error)),
    n_heldout = nrow(d)
  )
}

split_keys <- c("analysis_dataset", "source_type", "metric", "interpretation", "train_n_max", "exponent", "train_r2")
metric_summary <- do.call(rbind, lapply(split(heldout, interaction(heldout[split_keys], drop = TRUE)), function(d) {
  cbind(d[1, split_keys], summary_one(d))
}))
rownames(metric_summary) <- NULL
metric_summary <- metric_summary[order(metric_summary$analysis_dataset, metric_summary$mean_abs_pct_error), ]

dataset_keys <- c("analysis_dataset", "source_type")
dataset_summary <- do.call(rbind, lapply(split(heldout, interaction(heldout[dataset_keys], drop = TRUE)), function(d) {
  cbind(d[1, dataset_keys], summary_one(d))
}))
rownames(dataset_summary) <- NULL
dataset_summary <- dataset_summary[order(dataset_summary$mean_abs_pct_error), ]

write.csv(pred, file.path(out_dir, "scaling_extrapolation_predictions.csv"), row.names = FALSE)
write.csv(heldout, file.path(out_dir, "scaling_extrapolation_heldout.csv"), row.names = FALSE)
write.csv(metric_summary, file.path(out_dir, "scaling_extrapolation_metric_summary.csv"), row.names = FALSE)
write.csv(dataset_summary, file.path(out_dir, "scaling_extrapolation_dataset_summary.csv"), row.names = FALSE)
write.csv(future_pred, file.path(out_dir, "scaling_extrapolation_future_predictions.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 12))

plot_df <- pred[pred$metric %in% c(
  "lambda1_over_edge",
  "top20_mean_over_edge",
  "sqrt_normalized_spike_excess",
  "n_spikes",
  "noise_edge",
  "raw_log_lambda1_over_median",
  "glm_eta_lambda1_over_raw_median"
), ]

short_metric_labels <- c(
  lambda1_over_edge = "lambda1 / edge",
  top20_mean_over_edge = "top20 mean / edge",
  sqrt_normalized_spike_excess = "sqrt spike excess",
  n_spikes = "MP-edge spikes",
  noise_edge = "MP upper edge",
  raw_log_lambda1_over_median = "raw lambda1 / median",
  glm_eta_lambda1_over_raw_median = "GLM eta lambda1 / raw median"
)
plot_df$metric_label <- factor(
  unname(short_metric_labels[plot_df$metric]),
  levels = unname(short_metric_labels[unique(plot_df$metric)])
)
plot_df$split <- factor(plot_df$split, levels = c("fit", "heldout"))

validation_plot <- ggplot(plot_df, aes(n_cells, observed)) +
  geom_ribbon(aes(ymin = pred_low80, ymax = pred_high80), fill = "#A6CEE3", alpha = 0.28) +
  geom_line(aes(y = predicted), color = "#1F78B4", linewidth = 0.75) +
  geom_point(aes(shape = split), size = 2.2, color = "#222222") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  facet_grid(metric_label ~ analysis_dataset, scales = "free_y") +
  labs(
    title = "Small-consortium scaling fits extrapolated to larger cell counts",
    subtitle = "Blue line is a log-log power-law fit on the smallest available n; points after the fit window are held out.",
    x = "Cells",
    y = "Observed metric",
    shape = "Split"
  )

ggsave(file.path(out_dir, "scaling_extrapolation_validation.png"), validation_plot, width = 13, height = 12, dpi = 220)

error_plot <- ggplot(heldout, aes(reorder(interpretation, ape), ape, fill = analysis_dataset)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Held-out extrapolation error by metric",
    x = NULL,
    y = "Absolute percent error",
    fill = "Dataset"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "scaling_extrapolation_heldout_error.png"), error_plot, width = 10, height = 7, dpi = 220)

print(dataset_summary)
print(metric_summary)
message("Wrote scaling extrapolation validation outputs to ", out_dir)
