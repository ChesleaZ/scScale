#!/usr/bin/env Rscript

local_lib <- normalizePath(".rlib/R", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
})

source("scripts/core/mp_bbp.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

data_path <- arg_value(
  "--data",
  "data/GSE123025/GSE123025_Single_myeloid_1922_cells_processed_data.csv.gz"
)
out_dir <- arg_value("--out-dir", "outputs/comparisons/gse123025_mp_bbp")
n_cells <- as.integer(arg_value("--n-cells", "600"))
n_genes <- as.integer(arg_value("--n-genes", "300"))
seed <- as.integer(arg_value("--seed", "20260611"))
mp_R <- as.integer(arg_value("--mp-R", "20"))
mp_p_sim <- as.integer(arg_value("--mp-p-sim", "500"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", data_path)
raw_counts <- read.csv(gzfile(data_path), check.names = FALSE)
gene_names <- raw_counts[[1]]
expr_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_names
counts0 <- as(expr_mat, "dgCMatrix")

set.seed(seed)
cell_totals <- Matrix::colSums(counts0)
eligible <- names(cell_totals)[cell_totals > 0]
if (length(eligible) < n_cells) stop("Not enough eligible cells.", call. = FALSE)
cell_use <- sample(eligible, n_cells)
counts <- counts0[, cell_use, drop = FALSE]

log_counts <- log1p(as.matrix(counts))
gene_var <- apply(log_counts, 1, stats::var)
hvg <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(n_genes, length(gene_var)))]
counts_hvg <- as.matrix(counts[hvg, , drop = FALSE])

build_matrix_pearson_counts <- function(counts_mat) {
  cell_totals <- colSums(counts_mat)
  gene_totals <- rowSums(counts_mat)
  overall_total <- sum(cell_totals)
  cell_totals[cell_totals == 0] <- 1
  gene_totals[gene_totals == 0] <- 1
  expected <- outer(gene_totals, cell_totals) / overall_total
  residual <- (counts_mat - expected) / sqrt(expected + 1e-4)
  residual[residual > 10] <- 10
  residual[residual < -10] <- -10
  residual
}

mp_median_sim_github <- function(c_ratio, R = 20, p_sim = 500, seed = 1) {
  set.seed(seed)
  n_sim <- max(round(p_sim / c_ratio), 20)
  meds <- replicate(R, {
    g <- matrix(rnorm(p_sim * n_sim), nrow = p_sim, ncol = n_sim)
    ev <- eigen((g %*% t(g)) / n_sim, only.values = TRUE)$values
    stats::median(ev)
  })
  stats::median(meds)
}

github_d_from_mu <- function(mu, J, n, p) {
  c_ratio <- p / n
  mu_tilde <- mu / J
  term <- mu_tilde - (1 + c_ratio)
  disc <- term^2 - 4 * c_ratio
  d <- rep(NA_real_, length(mu))
  ok <- is.finite(disc) & disc >= 0
  d[ok] <- sqrt(pmax((term[ok] + sqrt(disc[ok])) / (2 * c_ratio), 0))
  d
}

x <- build_matrix_pearson_counts(counts_hvg)
n_real <- ncol(x)
p_real <- nrow(x)
c_ratio <- p_real / n_real

sv <- svd(t(x), nu = 0, nv = 0)$d
mu <- sort((sv^2) / n_real, decreasing = TRUE)

message("Matrix: p=", p_real, ", n=", n_real, ", c=", signif(c_ratio, 4))

mp_median <- mp_median_sim_github(c_ratio, R = mp_R, p_sim = mp_p_sim, seed = seed + 1)
J_github <- stats::median(mu) / mp_median
mp_upper_github <- J_github * (1 + sqrt(c_ratio))^2
d_github <- github_d_from_mu(mu, J_github, n_real, p_real)
d_github[mu <= mp_upper_github] <- NA_real_

fit_fixed <- fit_mp_robust(
  eigenvalues = mu,
  n_cells = n_real,
  n_features = p_real,
  fit_type = "fixed_gamma",
  edge_buffer = 1
)
fit_free <- fit_mp_robust(
  eigenvalues = mu,
  n_cells = n_real,
  n_features = p_real,
  fit_type = "free_gamma",
  edge_buffer = 1
)

spikes_fixed <- fit_bbp_spikes(
  eigenvalues = mu,
  n_cells = n_real,
  n_features = p_real,
  fit_type = "fixed_gamma",
  edge_buffer = 1
)

our_spikes <- spikes_fixed$spikes
our_spikes$github_d2_equivalent <- (our_spikes$population_eigenvalue_hat / fit_fixed$sigma2 - 1) / c_ratio

top_k <- seq_len(min(20, length(mu)))
top_df <- data.frame(
  rank = top_k,
  eigenvalue = mu[top_k],
  github_is_spike = mu[top_k] > mp_upper_github,
  github_d2 = d_github[top_k]^2,
  our_is_spike_fixed = mu[top_k] > fit_fixed$mp_upper,
  our_population_eigenvalue_fixed = NA_real_,
  our_d2_equivalent_fixed = NA_real_
)
if (nrow(our_spikes)) {
  hit <- match(top_df$eigenvalue, our_spikes$sample_eigenvalue)
  ok <- !is.na(hit)
  top_df$our_population_eigenvalue_fixed[ok] <- our_spikes$population_eigenvalue_hat[hit[ok]]
  top_df$our_d2_equivalent_fixed[ok] <- our_spikes$github_d2_equivalent[hit[ok]]
}

summary_df <- data.frame(
  method = c("github_mp_median", "ours_fixed_gamma_quantile", "ours_free_gamma_quantile"),
  gamma = c(c_ratio, fit_fixed$gamma, fit_free$gamma),
  sigma2_or_J = c(J_github, fit_fixed$sigma2, fit_free$sigma2),
  mp_upper = c(mp_upper_github, fit_fixed$mp_upper, fit_free$mp_upper),
  threshold = c(mp_upper_github, fit_fixed$mp_threshold, fit_free$mp_threshold),
  n_spikes = c(sum(mu > mp_upper_github), fit_fixed$n_outliers_above_mp, fit_free$n_outliers_above_mp),
  lambda1_over_edge = c(mu[1] / mp_upper_github, fit_fixed$lambda1_over_mp_upper, fit_free$lambda1_over_mp_upper),
  objective = c(NA_real_, fit_fixed$objective, fit_free$objective)
)

write.csv(summary_df, file.path(out_dir, "mp_fit_summary.csv"), row.names = FALSE)
write.csv(top_df, file.path(out_dir, "top_eigen_spike_comparison.csv"), row.names = FALSE)
write.csv(data.frame(rank = seq_along(mu), eigenvalue = mu), file.path(out_dir, "eigenvalues.csv"), row.names = FALSE)

png(file.path(out_dir, "mp_fit_comparison.png"), width = 1400, height = 900, res = 180)
hist(mu, breaks = 80, col = "grey85", border = "white", main = "GSE123025 Pearson residual covariance spectrum", xlab = "Eigenvalue")
abline(v = mp_upper_github, col = "#D55E00", lwd = 2)
abline(v = fit_fixed$mp_upper, col = "#0072B2", lwd = 2)
abline(v = fit_free$mp_upper, col = "#009E73", lwd = 2)
legend(
  "topright",
  legend = c("GitHub MP median edge", "Ours fixed-gamma edge", "Ours free-gamma edge"),
  col = c("#D55E00", "#0072B2", "#009E73"),
  lwd = 2,
  bty = "n"
)
dev.off()

html <- c(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>MP/BBP comparison</title>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:28px;max-width:1100px} table{border-collapse:collapse} th,td{border:1px solid #ddd;padding:6px 8px;text-align:right} th{text-align:center;background:#f4f6f8} td:first-child{text-align:left} img{max-width:100%;border:1px solid #ddd}</style>",
  "</head><body>",
  "<h1>GSE123025 MP/BBP comparison</h1>",
  sprintf("<p>Same Pearson residual matrix: %d genes x %d cells. HVGs selected by naive log1p variance.</p>", p_real, n_real),
  "<h2>Fit Summary</h2>",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "<br>"),
  "<h2>Top Eigenvalues</h2>",
  paste(capture.output(print(top_df, row.names = FALSE)), collapse = "<br>"),
  "<h2>Spectrum</h2><img src=\"mp_fit_comparison.png\">",
  "</body></html>"
)
writeLines(html, file.path(out_dir, "index.html"))

message("Wrote comparison to ", out_dir)
print(summary_df)
print(top_df)
