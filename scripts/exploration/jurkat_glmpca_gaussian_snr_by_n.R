#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: jurkat_glmpca_gaussian_snr_by_n.R <counts_csv> <out_dir> [dataset_name]", call. = FALSE)
}

counts_csv <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "jurkat")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(getwd(), ".rlib"), .libPaths()))

library(glmpca)
library(ggplot2)

set.seed(as.integer(Sys.getenv("SNR_SEED", "1")))

n_factors_requested <- as.integer(Sys.getenv("SNR_FACTORS", "20"))
max_iter <- as.integer(Sys.getenv("SNR_GLMPCA_MAX_ITER", "160"))
min_iter <- min(20, max_iter - 1)
cell_counts_requested <- as.integer(strsplit(Sys.getenv("SNR_CELL_COUNTS", "250,500,1000,2000,5000"), ",")[[1]])

message("Reading sampled count block from ", counts_csv)
counts_df <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
Y_full <- as.matrix(counts_df)
storage.mode(Y_full) <- "integer"

trace_cov <- function(X) {
  X <- as.matrix(X)
  Xc <- sweep(X, 1, rowMeans(X), "-")
  sv <- svd(Xc, nu = 0, nv = 0)$d
  sum(sv^2) / (ncol(Xc) - 1)
}

log_normalized <- function(Y) {
  lib_size <- colSums(Y)
  lib_size[lib_size <= 0] <- 1
  log1p(sweep(Y, 2, lib_size, "/") * 1e4)
}

gaussian_snr <- function(Y, rank) {
  X <- log_normalized(Y)
  Xc <- sweep(X, 1, rowMeans(X), "-")
  eig <- svd(Xc, nu = 0, nv = 0)$d^2 / (ncol(Xc) - 1)
  rank <- min(rank, length(eig) - 1)
  signal <- sum(eig[seq_len(rank)])
  noise <- sum(eig[(rank + 1):length(eig)])
  data.frame(
    method = "gaussian_log_pca",
    scale = "log_normalized",
    signal_trace = signal,
    raw_trace = signal + noise,
    noise_trace = noise,
    snr = signal / max(noise, 1e-12),
    signal_fraction = signal / max(signal + noise, 1e-12),
    lambda1 = eig[1],
    top_rank_sum = signal,
    glmpca_final_deviance = NA_real_,
    glmpca_final_deviance_delta = NA_real_,
    glmpca_relative_deviance_delta = NA_real_
  )
}

reconstruct_glmpca_mean <- function(fit) {
  factors <- as.matrix(fit$factors)
  loadings <- as.matrix(fit$loadings)
  X <- as.matrix(fit$X)
  coefX <- as.matrix(fit$coefX)

  eta <- loadings %*% t(factors)
  if (!is.null(X) && !is.null(coefX)) {
    eta <- eta + coefX %*% t(X)
  }
  if (!is.null(fit$offsets)) {
    eta <- sweep(eta, 2, fit$offsets, "+")
  }

  exp(pmin(eta, 700))
}

glmpca_snr <- function(Y, rank) {
  fit <- glmpca(
    Y,
    L = rank,
    fam = "poi",
    ctl = list(minIter = min_iter, maxIter = max_iter, verbose = FALSE)
  )

  M_hat <- reconstruct_glmpca_mean(fit)
  signal <- trace_cov(M_hat)
  raw <- trace_cov(Y)
  noise <- raw - signal

  data.frame(
    method = "poisson_glmpca_mean",
    scale = "counts",
    signal_trace = signal,
    raw_trace = raw,
    noise_trace = noise,
    snr = signal / max(noise, 1e-12),
    signal_fraction = signal / max(raw, 1e-12),
    lambda1 = NA_real_,
    top_rank_sum = NA_real_,
    glmpca_final_deviance = tail(fit$dev, 1),
    glmpca_final_deviance_delta = ifelse(length(fit$dev) > 1, fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1), NA_real_),
    glmpca_relative_deviance_delta = ifelse(length(fit$dev) > 1, (fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1)) / abs(tail(fit$dev, 1)), NA_real_)
  )
}

n_total <- ncol(Y_full)
cell_counts <- sort(unique(pmin(cell_counts_requested, n_total)))
cell_counts <- cell_counts[cell_counts > n_factors_requested + 2]

cell_order <- sample.int(n_total, n_total, replace = FALSE)

rows <- list()
for (n_cells in cell_counts) {
  message("Analyzing n_cells = ", n_cells)
  cols <- sort(cell_order[seq_len(n_cells)])
  Y <- Y_full[, cols, drop = FALSE]
  rank <- min(n_factors_requested, nrow(Y) - 1, ncol(Y) - 2)

  g_snr <- gaussian_snr(Y, rank)
  p_snr <- glmpca_snr(Y, rank)

  rows[[as.character(n_cells)]] <- rbind(
    cbind(
      dataset = dataset_name,
      n_cells = n_cells,
      n_genes = nrow(Y),
      rank = rank,
      mean_umi_per_cell_hvg = mean(colSums(Y)),
      g_snr
    ),
    cbind(
      dataset = dataset_name,
      n_cells = n_cells,
      n_genes = nrow(Y),
      rank = rank,
      mean_umi_per_cell_hvg = mean(colSums(Y)),
      p_snr
    )
  )
}

snr_df <- do.call(rbind, rows)
rownames(snr_df) <- NULL

numeric_cols <- setdiff(colnames(snr_df), c("dataset", "method", "scale"))
for (col in numeric_cols) {
  snr_df[[col]] <- as.numeric(snr_df[[col]])
}

out_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_gaussian_snr_by_n.csv"))
write.csv(snr_df, out_file, row.names = FALSE)

snr_df$method_label <- factor(
  snr_df$method,
  levels = c("gaussian_log_pca", "poisson_glmpca_mean"),
  labels = c("Gaussian log-PCA", "Poisson GLM-PCA mean")
)

theme_set(theme_classic(base_size = 13))

snr_plot <- ggplot(snr_df, aes(n_cells, snr, color = method_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": empirical SNR scaling with cells"),
    x = "Cells",
    y = "SNR (log scale)",
    color = "Method"
  )

fraction_plot <- ggplot(snr_df, aes(n_cells, signal_fraction, color = method_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = paste0(dataset_name, ": signal fraction scaling with cells"),
    x = "Cells",
    y = "Signal trace / raw trace",
    color = "Method"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_glmpca_gaussian_snr_by_n.png")), snr_plot, width = 8, height = 5, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_glmpca_gaussian_signal_fraction_by_n.png")), fraction_plot, width = 8, height = 5, dpi = 220)

print(snr_df)
message("Wrote SNR table to ", out_file)
