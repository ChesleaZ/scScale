#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: glmpca_eigenvalue_comparison.R <counts_csv> <out_dir> [dataset_name]", call. = FALSE)
}

counts_csv <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "jurkat")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(getwd(), ".rlib"), .libPaths()))

library(glmpca)
library(ggplot2)
library(Matrix)

set.seed(1)

n_factors <- as.integer(Sys.getenv("GLMPCA_FACTORS", "25"))
max_iter <- as.integer(Sys.getenv("GLMPCA_MAX_ITER", "80"))
min_iter <- min(20, max_iter - 1)
depths <- as.numeric(strsplit(Sys.getenv("GLMPCA_DEPTHS", "0.1,0.25,0.5,1"), ",")[[1]])

message("Reading sampled count block from ", counts_csv)
counts_df <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
Y_full <- as.matrix(counts_df)
storage.mode(Y_full) <- "integer"

thin_counts <- function(Y, rate) {
  if (rate >= 1) {
    return(Y)
  }
  matrix(rbinom(length(Y), size = as.vector(Y), prob = rate), nrow = nrow(Y), ncol = ncol(Y))
}

log_pca_eigenvalues <- function(Y) {
  lib_size <- colSums(Y)
  lib_size[lib_size <= 0] <- 1
  X <- log1p(sweep(Y, 2, lib_size, "/") * 1e4)
  X <- sweep(X, 1, rowMeans(X), "-")
  sv <- svd(X, nu = 0, nv = 0)$d
  eig <- sv^2 / (ncol(X) - 1)
  if (length(eig) < nrow(X)) {
    eig <- c(eig, rep(0, nrow(X) - length(eig)))
  }
  sort(eig, decreasing = TRUE)
}

glmpca_eta_eigenvalues <- function(fit) {
  F <- as.matrix(fit$factors)
  L <- as.matrix(fit$loadings)
  eta <- t(L %*% t(F))
  eta <- sweep(eta, 2, colMeans(eta), "-")
  sv <- svd(eta, nu = 0, nv = 0)$d
  sv^2 / (nrow(eta) - 1)
}

all_eigs <- list()
all_summary <- list()

for (rate in depths) {
  message("Fitting depth rate = ", rate)
  Y <- thin_counts(Y_full, rate)
  mean_umi <- mean(colSums(Y))

  pca_eig <- log_pca_eigenvalues(Y)
  pca_eig <- pca_eig[seq_len(min(length(pca_eig), n_factors))]

  fit <- glmpca(
    Y,
    L = n_factors,
    fam = "poi",
    ctl = list(minIter = min_iter, maxIter = max_iter, verbose = FALSE)
  )
  glm_eig <- glmpca_eta_eigenvalues(fit)
  glm_eig <- glm_eig[seq_len(min(length(glm_eig), n_factors))]

  all_eigs[[as.character(rate)]] <- rbind(
    data.frame(
      dataset = dataset_name,
      depth_rate = rate,
      mean_umi_per_cell = mean_umi,
      method = "log_normalized_pca",
      rank = seq_along(pca_eig),
      eigenvalue = pca_eig
    ),
    data.frame(
      dataset = dataset_name,
      depth_rate = rate,
      mean_umi_per_cell = mean_umi,
      method = "poisson_glmpca_eta",
      rank = seq_along(glm_eig),
      eigenvalue = glm_eig
    )
  )

  all_summary[[as.character(rate)]] <- data.frame(
    dataset = dataset_name,
    depth_rate = rate,
    mean_umi_per_cell = mean_umi,
    n_genes = nrow(Y),
    n_cells = ncol(Y),
    n_factors = n_factors,
    max_iter = max_iter,
    glmpca_final_deviance = tail(fit$dev, 1),
    glmpca_final_deviance_delta = ifelse(length(fit$dev) > 1, fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1), NA_real_),
    glmpca_relative_deviance_delta = ifelse(length(fit$dev) > 1, (fit$dev[length(fit$dev) - 1] - tail(fit$dev, 1)) / abs(tail(fit$dev, 1)), NA_real_),
    log_pca_lambda1 = pca_eig[1],
    glmpca_eta_lambda1 = glm_eig[1],
    log_pca_top10_sum = sum(head(pca_eig, 10)),
    glmpca_eta_top10_sum = sum(head(glm_eig, 10))
  )
}

eig_df <- do.call(rbind, all_eigs)
summary_df <- do.call(rbind, all_summary)
rownames(eig_df) <- NULL
rownames(summary_df) <- NULL

eig_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_eigenvalue_comparison.csv"))
summary_file <- file.path(out_dir, paste0(dataset_name, "_glmpca_eigenvalue_summary.csv"))
write.csv(eig_df, eig_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)

eig_df$method_label <- factor(
  eig_df$method,
  levels = c("log_normalized_pca", "poisson_glmpca_eta"),
  labels = c("log-normalized PCA", "Poisson GLM-PCA eta")
)
eig_df$depth_label <- paste0(
  "rate ", eig_df$depth_rate,
  " (mean UMI ", format(round(eig_df$mean_umi_per_cell), big.mark = ","), ")"
)
eig_df$depth_label <- factor(eig_df$depth_label, levels = unique(eig_df$depth_label))

theme_set(theme_classic(base_size = 13))

scree_plot <- ggplot(eig_df, aes(rank, eigenvalue, color = method_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_y_log10() +
  facet_wrap(~ depth_label, ncol = 2) +
  labs(
    title = paste0(dataset_name, ": GLM-PCA vs log-PCA eigenvalues across count depth"),
    x = "Rank",
    y = "Eigenvalue (log scale)",
    color = "Method"
  )

top_plot <- ggplot(summary_df, aes(mean_umi_per_cell)) +
  geom_line(aes(y = log_pca_lambda1, color = "log-normalized PCA"), linewidth = 0.9) +
  geom_point(aes(y = log_pca_lambda1, color = "log-normalized PCA"), size = 2) +
  geom_line(aes(y = glmpca_eta_lambda1, color = "Poisson GLM-PCA eta"), linewidth = 0.9) +
  geom_point(aes(y = glmpca_eta_lambda1, color = "Poisson GLM-PCA eta"), size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": top eigenvalue vs mean UMI per cell"),
    x = "Mean UMI per cell",
    y = "Top eigenvalue (log scale)",
    color = "Method"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_glmpca_scree_by_depth.png")), scree_plot, width = 10, height = 7, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_glmpca_top_eigenvalue_by_depth.png")), top_plot, width = 8, height = 5, dpi = 220)

print(summary_df)
message("Wrote eigenvalues to ", eig_file)
message("Wrote summary to ", summary_file)
