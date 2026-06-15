#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: save_glmpca_log_eigenvalues_by_n.R <counts_csv> <out_dir> [dataset_name]", call. = FALSE)
}

counts_csv <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(getwd(), ".rlib"), .libPaths()))

library(glmpca)

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
  as.matrix(fit$loadings) %*% t(as.matrix(fit$factors))
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

n_total <- ncol(Y_full)
cell_counts <- sort(unique(pmin(cell_counts_requested, n_total)))
cell_counts <- cell_counts[cell_counts > n_factors_requested + 2]
cell_order <- sample.int(n_total, n_total, replace = FALSE)

eig_rows <- list()
summary_rows <- list()

for (n_cells in cell_counts) {
  message("Analyzing n_cells = ", n_cells)
  cols <- sort(cell_order[seq_len(n_cells)])
  Y <- Y_full[, cols, drop = FALSE]
  rank <- min(n_factors_requested, nrow(Y) - 1, ncol(Y) - 2)

  X_log <- log_normalized(Y)
  raw_log_eig <- cov_eigenvalues(X_log)

  init <- if (init_method == "svd") svd_init(Y, rank) else list(factors = NULL, loadings = NULL)
  fit <- glmpca(
    Y,
    L = rank,
    fam = "poi",
    init = init,
    ctl = list(minIter = min_iter, maxIter = max_iter, verbose = FALSE)
  )
  eta_eig <- cov_eigenvalues(eta_lowrank(fit))
  eta_eig <- eta_eig[seq_len(min(length(eta_eig), rank))]

  eig_rows[[length(eig_rows) + 1]] <- data.frame(
    dataset = dataset_name,
    n_cells = n_cells,
    n_genes = nrow(Y),
    rank_fit = rank,
    method = "raw_log_normalized",
    eig_rank = seq_along(raw_log_eig),
    eigenvalue = raw_log_eig,
    log10_eigenvalue = log10(pmax(raw_log_eig, .Machine$double.xmin))
  )

  eig_rows[[length(eig_rows) + 1]] <- data.frame(
    dataset = dataset_name,
    n_cells = n_cells,
    n_genes = nrow(Y),
    rank_fit = rank,
    method = "glmpca_eta_lowrank",
    eig_rank = seq_along(eta_eig),
    eigenvalue = eta_eig,
    log10_eigenvalue = log10(pmax(eta_eig, .Machine$double.xmin))
  )

  raw_median <- median(raw_log_eig)
  raw_tail_median <- median(raw_log_eig[(rank + 1):length(raw_log_eig)])
  eta_positive <- eta_eig[eta_eig > 1e-10]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    dataset = dataset_name,
    n_cells = n_cells,
    n_genes = nrow(Y),
    rank_fit = rank,
    mean_umi_per_cell_hvg = mean(colSums(Y)),
    glmpca_init = init_method,
    glmpca_final_deviance = tail(fit$dev, 1),
    glmpca_final_deviance_delta = ifelse(length(fit$dev) > 1, fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1), NA_real_),
    glmpca_relative_deviance_delta = ifelse(length(fit$dev) > 1, (fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1)) / abs(tail(fit$dev, 1)), NA_real_),
    raw_log_lambda1 = raw_log_eig[1],
    raw_log_median = raw_median,
    raw_log_tail_median = raw_tail_median,
    raw_log_lambda1_over_median = raw_log_eig[1] / raw_median,
    raw_log_lambda1_over_tail_median = raw_log_eig[1] / raw_tail_median,
    glm_eta_lambda1 = eta_eig[1],
    glm_eta_nonzero_median = median(eta_positive),
    glm_eta_lambda1_over_raw_median = eta_eig[1] / raw_median,
    glm_eta_lambda1_over_raw_tail_median = eta_eig[1] / raw_tail_median,
    glm_eta_lambda1_over_eta_nonzero_median = eta_eig[1] / median(eta_positive)
  )
}

eig_df <- do.call(rbind, eig_rows)
summary_df <- do.call(rbind, summary_rows)
rownames(eig_df) <- NULL
rownames(summary_df) <- NULL

eig_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_log_eigenvalues_by_n.csv"))
summary_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_log_eigenvalues_summary_by_n.csv"))

write.csv(eig_df, eig_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)

print(summary_df)
message("Wrote reusable eigenvalues to ", eig_file)
message("Wrote summary to ", summary_file)
