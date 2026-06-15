#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: glmpca_spike_bulk_proxy_by_n.R <counts_csv> <out_dir> [dataset_name]", call. = FALSE)
}

counts_csv <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(getwd(), ".rlib"), .libPaths()))

library(glmpca)
library(ggplot2)

set.seed(as.integer(Sys.getenv("SNR_SEED", "1")))

n_factors_requested <- as.integer(Sys.getenv("SNR_FACTORS", "20"))
max_iter <- as.integer(Sys.getenv("SNR_GLMPCA_MAX_ITER", "220"))
min_iter <- min(20, max_iter - 1)
cell_counts_requested <- as.integer(strsplit(Sys.getenv("SNR_CELL_COUNTS", "250,500,1000,2000,5000"), ",")[[1]])
init_method <- Sys.getenv("SNR_GLMPCA_INIT", "svd")

message("Reading sampled count block from ", counts_csv)
counts_df <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
Y_full <- as.matrix(counts_df)
storage.mode(Y_full) <- "integer"

cov_eigenvalues <- function(X) {
  X <- as.matrix(X)
  Xc <- sweep(X, 1, rowMeans(X), "-")
  eig <- svd(Xc, nu = 0, nv = 0)$d^2 / (ncol(Xc) - 1)
  sort(eig, decreasing = TRUE)
}

log_normalized <- function(Y) {
  lib_size <- colSums(Y)
  lib_size[lib_size <= 0] <- 1
  log1p(sweep(Y, 2, lib_size, "/") * 1e4)
}

eta_lowrank <- function(fit) {
  loadings <- as.matrix(fit$loadings)
  factors <- as.matrix(fit$factors)
  loadings %*% t(factors)
}

svd_init <- function(Y, rank) {
  X <- log_normalized(Y)
  X <- sweep(X, 1, rowMeans(X), "-")
  sv <- svd(X, nu = rank, nv = rank)
  d_sqrt <- sqrt(pmax(sv$d[seq_len(rank)], 0))
  list(
    loadings = sweep(sv$u[, seq_len(rank), drop = FALSE], 2, d_sqrt, "*"),
    factors = sweep(sv$v[, seq_len(rank), drop = FALSE], 2, d_sqrt, "*")
  )
}

analyze_subset <- function(Y, rank) {
  raw_log_eig <- cov_eigenvalues(log_normalized(Y))
  raw_log_bulk_median <- median(raw_log_eig)
  raw_log_tail_median <- median(raw_log_eig[(rank + 1):length(raw_log_eig)])

  init <- if (init_method == "svd") svd_init(Y, rank) else list(factors = NULL, loadings = NULL)
  fit <- glmpca(
    Y,
    L = rank,
    fam = "poi",
    init = init,
    ctl = list(minIter = min_iter, maxIter = max_iter, verbose = FALSE)
  )

  eta_eig <- cov_eigenvalues(eta_lowrank(fit))
  eta_positive <- eta_eig[eta_eig > 1e-10]

  data.frame(
    raw_log_lambda1 = raw_log_eig[1],
    raw_log_top5_mean = mean(head(raw_log_eig, 5)),
    raw_log_bulk_median = raw_log_bulk_median,
    raw_log_tail_median = raw_log_tail_median,
    glmpca_init = init_method,
    raw_log_lambda1_over_median = raw_log_eig[1] / raw_log_bulk_median,
    raw_log_lambda1_over_tail_median = raw_log_eig[1] / raw_log_tail_median,
    glm_eta_lambda1 = eta_eig[1],
    glm_eta_top5_mean = mean(head(eta_eig, min(5, length(eta_eig)))),
    glm_eta_nonzero_median = median(eta_positive),
    glm_eta_lambda1_over_raw_median = eta_eig[1] / raw_log_bulk_median,
    glm_eta_lambda1_over_raw_tail_median = eta_eig[1] / raw_log_tail_median,
    glm_eta_lambda1_over_eta_nonzero_median = eta_eig[1] / median(eta_positive),
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
  rows[[as.character(n_cells)]] <- cbind(
    dataset = dataset_name,
    n_cells = n_cells,
    n_genes = nrow(Y),
    rank = rank,
    mean_umi_per_cell_hvg = mean(colSums(Y)),
    analyze_subset(Y, rank)
  )
}

proxy_df <- do.call(rbind, rows)
rownames(proxy_df) <- NULL
numeric_cols <- setdiff(colnames(proxy_df), c("dataset", "glmpca_init"))
for (col in numeric_cols) {
  proxy_df[[col]] <- as.numeric(proxy_df[[col]])
}

out_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_spike_bulk_proxy_by_n.csv"))
write.csv(proxy_df, out_file, row.names = FALSE)

plot_df <- rbind(
  data.frame(
    n_cells = proxy_df$n_cells,
    proxy = proxy_df$raw_log_lambda1_over_median,
    proxy_type = "Raw log-PCA lambda1 / raw median"
  ),
  data.frame(
    n_cells = proxy_df$n_cells,
    proxy = proxy_df$glm_eta_lambda1_over_raw_median,
    proxy_type = "GLM eta lambda1 / raw median"
  ),
  data.frame(
    n_cells = proxy_df$n_cells,
    proxy = proxy_df$glm_eta_lambda1_over_eta_nonzero_median,
    proxy_type = "GLM eta lambda1 / eta nonzero median"
  )
)

plot_df$proxy_type <- factor(
  plot_df$proxy_type,
  levels = c(
    "Raw log-PCA lambda1 / raw median",
    "GLM eta lambda1 / raw median",
    "GLM eta lambda1 / eta nonzero median"
  )
)

theme_set(theme_classic(base_size = 13))

proxy_plot <- ggplot(plot_df, aes(n_cells, proxy, color = proxy_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": spike-to-bulk proxy scaling"),
    x = "Cells",
    y = "Spike / bulk proxy (log scale)",
    color = "Proxy"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_glmpca_spike_bulk_proxy_by_n.png")), proxy_plot, width = 9, height = 5.5, dpi = 220)

print(proxy_df)
message("Wrote spike-to-bulk proxy table to ", out_file)
