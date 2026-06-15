#!/usr/bin/env Rscript

local_lib <- normalizePath(".rlib/R", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
})

source("scripts/core/mp_bbp.R")
source("scripts/core/spectral_overlap_mi.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

data_path <- arg_value("--data", "data/GSE123025/GSE123025_Single_myeloid_1922_cells_processed_data.csv.gz")
out_dir <- arg_value("--out-dir", "outputs/comparisons/gse123025_mi")
n_cells <- as.integer(arg_value("--n-cells", "400"))
n_genes <- as.integer(arg_value("--n-genes", "300"))
seed <- as.integer(arg_value("--seed", "20260611"))
u_ref <- as.numeric(arg_value("--u-ref", "300000"))
u_grid <- as.numeric(strsplit(arg_value("--u-grid", "1000,10000,50000,100000,150000"), ",", fixed = TRUE)[[1]])
r_mi <- as.integer(arg_value("--r-mi", "5"))
mp_R <- as.integer(arg_value("--mp-R", "10"))
mp_p_sim <- as.integer(arg_value("--mp-p-sim", "300"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

downsample_counts_to_target <- function(counts, target_umi, seed = 1) {
  set.seed(seed)
  counts <- as(counts, "dgCMatrix")
  col_sums <- Matrix::colSums(counts)
  to_down <- which(col_sums > target_umi)
  if (length(to_down)) {
    probs <- target_umi / col_sums[to_down]
    names(probs) <- as.character(to_down)
    for (j in to_down) {
      start <- counts@p[j] + 1
      end <- counts@p[j + 1]
      if (start <= end) {
        counts@x[start:end] <- stats::rbinom(
          n = end - start + 1,
          size = counts@x[start:end],
          prob = probs[as.character(j)]
        )
      }
    }
  }
  drop0(counts)
}

sparse_nonneg_subtract <- function(a, b) {
  c <- as(a, "dgCMatrix") - as(b, "dgCMatrix")
  c@x[c@x < 0] <- 0
  drop0(c)
}

build_matrix_pearson_counts <- function(counts_mat) {
  counts_mat <- as.matrix(counts_mat)
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

mp_median_sim_github <- function(c_ratio, R = 10, p_sim = 300, seed = 1) {
  set.seed(seed)
  n_sim <- max(round(p_sim / c_ratio), 20)
  meds <- replicate(R, {
    g <- matrix(rnorm(p_sim * n_sim), nrow = p_sim, ncol = n_sim)
    ev <- eigen((g %*% t(g)) / n_sim, only.values = TRUE)$values
    stats::median(ev)
  })
  stats::median(meds)
}

github_d_from_mu <- function(mu, j_hat, n, p) {
  c_ratio <- p / n
  mu_tilde <- mu / j_hat
  term <- mu_tilde - (1 + c_ratio)
  disc <- term^2 - 4 * c_ratio
  d <- rep(NA_real_, length(mu))
  ok <- is.finite(disc) & disc >= 0
  d[ok] <- sqrt(pmax((term[ok] + sqrt(disc[ok])) / (2 * c_ratio), 0))
  d
}

estimate_github_spikes <- function(x, j_hat, r_keep = r_mi) {
  n <- ncol(x)
  p <- nrow(x)
  c_ratio <- p / n
  sv <- svd(t(x), nu = 0, nv = 0)$d
  mu <- sort((sv^2) / n, decreasing = TRUE)
  edge <- j_hat * (1 + sqrt(c_ratio))^2
  d <- github_d_from_mu(mu, j_hat, n, p)
  d[mu <= edge] <- NA_real_
  list(mu = mu, edge = edge, d = d, n_spikes = sum(mu > edge, na.rm = TRUE))
}

estimate_our_spikes <- function(x, fit_type = c("median_fixed_gamma", "fixed_gamma", "free_gamma")) {
  fit_type <- match.arg(fit_type)
  n <- ncol(x)
  p <- nrow(x)
  sv <- svd(t(x), nu = 0, nv = 0)$d
  mu <- sort((sv^2) / n, decreasing = TRUE)
  fit <- fit_bbp_spikes(mu, n_cells = n, n_features = p, fit_type = fit_type, edge_buffer = 1)
  d_equiv <- rep(NA_real_, length(mu))
  if (nrow(fit$spikes)) {
    hit <- match(fit$spikes$sample_eigenvalue, mu)
    d_equiv[hit] <- sqrt(pmax((fit$spikes$population_eigenvalue_hat / fit$bulk_fit$sigma2 - 1) / (p / n), 0))
  }
  list(mu = mu, edge = fit$bulk_fit$mp_upper, d = d_equiv, n_spikes = nrow(fit$spikes), fit = fit$bulk_fit)
}

github_empirical_overlap_ref_cells <- function(x_cur, x_ref, r = r_mi) {
  common_genes <- intersect(rownames(x_cur), rownames(x_ref))
  common_cells <- intersect(colnames(x_cur), colnames(x_ref))
  if (length(common_genes) < 2 || length(common_cells) < 2) {
    return(list(mi = NA_real_, gamma = rep(NA_real_, r), r_eff = 0))
  }
  xc <- x_cur[common_genes, common_cells, drop = FALSE]
  xr <- x_ref[common_genes, common_cells, drop = FALSE]
  r_eff <- min(r, ncol(xc), nrow(xc), ncol(xr), nrow(xr))
  v_hat <- svd(xc, nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  v_ref <- svd(xr, nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  s <- t(v_hat) %*% v_ref
  gamma <- svd(s, nu = 0, nv = 0)$d^2
  gamma <- pmin(pmax(gamma, 0), 1 - 1e-12)
  out <- rep(NA_real_, r)
  out[seq_along(gamma)] <- gamma
  list(mi = -0.5 * sum(log(1 - gamma)), gamma = out, r_eff = r_eff)
}

github_theoretical_overlap_ref_cells <- function(d_cur, d_ref, c_cur, c_ref = c_cur, r = r_mi) {
  gamma_one <- function(d, c_ratio) {
    d_use <- d[is.finite(d) & d > 0]
    d_use <- d_use[seq_len(min(r, length(d_use)))]
    out <- rep(NA_real_, r)
    if (!length(d_use)) return(out)
    g <- (d_use^4 * c_ratio - 1) / (d_use^2 * (c_ratio * d_use^2 + 1))
    g[d_use^4 <= 1 / c_ratio] <- 0
    out[seq_along(g)] <- pmin(pmax(g, 0), 1 - 1e-12)
    out
  }
  gamma_cur <- gamma_one(d_cur, c_cur)
  gamma_ref <- gamma_one(d_ref, c_ref)
  gamma_cur[!is.finite(gamma_cur)] <- 0
  gamma_ref[!is.finite(gamma_ref)] <- 0
  gamma_double <- pmin(pmax(gamma_cur * gamma_ref, 0), 1 - 1e-12)
  gamma_use <- gamma_double[is.finite(gamma_double) & gamma_double > 0]
  list(
    mi = if (!length(gamma_use)) 0 else -0.5 * sum(log(1 - gamma_use)),
    gamma_current = gamma_cur,
    gamma_reference = gamma_ref,
    gamma_double = gamma_double
  )
}

message("Reading ", data_path)
raw_counts <- read.csv(gzfile(data_path), check.names = FALSE)
gene_names <- raw_counts[[1]]
expr_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_names
counts0 <- as(expr_mat, "dgCMatrix")

set.seed(seed)
min_total_needed <- u_ref + max(u_grid) + 2000
cell_totals <- Matrix::colSums(counts0)
eligible <- names(cell_totals)[cell_totals >= min_total_needed]
if (length(eligible) < n_cells) {
  stop("Need ", n_cells, " cells with at least ", min_total_needed, " UMI; found ", length(eligible), call. = FALSE)
}
cell_use <- sample(eligible, n_cells)
counts_full <- counts0[, cell_use, drop = FALSE]

log_counts <- log1p(as.matrix(counts_full))
gene_var <- apply(log_counts, 1, stats::var)
genes_use <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(n_genes, length(gene_var)))]

counts_ref <- downsample_counts_to_target(counts_full, u_ref, seed = seed + 10)
counts_pool <- sparse_nonneg_subtract(counts_full, counts_ref)
x_ref <- build_matrix_pearson_counts(counts_ref[genes_use, , drop = FALSE])

n_ref <- ncol(x_ref)
p_ref <- nrow(x_ref)
c_ref <- p_ref / n_ref
sv_ref <- svd(t(x_ref), nu = 0, nv = 0)$d
mu_ref <- sort((sv_ref^2) / n_ref, decreasing = TRUE)
j_ref <- stats::median(mu_ref) / mp_median_sim_github(c_ref, R = mp_R, p_sim = mp_p_sim, seed = seed + 20)

ref_github <- estimate_github_spikes(x_ref, j_ref)
ref_ours_median <- estimate_our_spikes(x_ref, "median_fixed_gamma")
ref_ours_fixed <- estimate_our_spikes(x_ref, "fixed_gamma")
ref_ours_free <- estimate_our_spikes(x_ref, "free_gamma")

rows <- lapply(seq_along(u_grid), function(i) {
  u <- u_grid[i]
  counts_cur <- downsample_counts_to_target(counts_pool, u, seed = seed + 100 + i)
  x_cur <- build_matrix_pearson_counts(counts_cur[genes_use, , drop = FALSE])
  n_cur <- ncol(x_cur)
  p_cur <- nrow(x_cur)
  c_cur <- p_cur / n_cur

  cur_github <- estimate_github_spikes(x_cur, j_ref)
  cur_ours_median <- estimate_our_spikes(x_cur, "median_fixed_gamma")
  cur_ours_fixed <- estimate_our_spikes(x_cur, "fixed_gamma")
  cur_ours_free <- estimate_our_spikes(x_cur, "free_gamma")

  emp_github <- github_empirical_overlap_ref_cells(x_cur, x_ref, r = r_mi)
  emp_ours <- cell_subspace_overlap_mi(x_cur, x_ref, r = r_mi)

  th_github <- github_theoretical_overlap_ref_cells(cur_github$d, ref_github$d, c_cur, c_ref, r = r_mi)
  th_ours_same_d <- theoretical_overlap_mi_cells(cur_github$d, ref_github$d, c_cur, c_ref, r = r_mi)
  th_ours_median <- theoretical_overlap_mi_cells(cur_ours_median$d, ref_ours_median$d, c_cur, c_ref, r = r_mi)
  th_ours_fixed <- theoretical_overlap_mi_cells(cur_ours_fixed$d, ref_ours_fixed$d, c_cur, c_ref, r = r_mi)
  th_ours_free <- theoretical_overlap_mi_cells(cur_ours_free$d, ref_ours_free$d, c_cur, c_ref, r = r_mi)

  data.frame(
    umi = u,
    n_cells = n_cur,
    n_genes = p_cur,
    c_ratio = c_cur,
    j_ref_github = j_ref,
    emp_mi_github = emp_github$mi,
    emp_mi_ours = emp_ours$mi,
    emp_mi_abs_diff = abs(emp_github$mi - emp_ours$mi),
    theory_mi_github = th_github$mi,
    theory_mi_ours_same_github_d = th_ours_same_d$mi,
    theory_same_d_abs_diff = abs(th_github$mi - th_ours_same_d$mi),
    theory_mi_ours_median_spikes = th_ours_median$mi,
    theory_mi_ours_fixed_spikes = th_ours_fixed$mi,
    theory_mi_ours_free_spikes = th_ours_free$mi,
    n_spikes_cur_github = cur_github$n_spikes,
    n_spikes_ref_github = ref_github$n_spikes,
    n_spikes_cur_ours_median = cur_ours_median$n_spikes,
    n_spikes_ref_ours_median = ref_ours_median$n_spikes,
    n_spikes_cur_ours_fixed = cur_ours_fixed$n_spikes,
    n_spikes_ref_ours_fixed = ref_ours_fixed$n_spikes,
    n_spikes_cur_ours_free = cur_ours_free$n_spikes,
    n_spikes_ref_ours_free = ref_ours_free$n_spikes
  )
})

result <- do.call(rbind, rows)
write.csv(result, file.path(out_dir, "mi_comparison_summary.csv"), row.names = FALSE)

png(file.path(out_dir, "mi_comparison_curves.png"), width = 1500, height = 950, res = 180)
ylim <- range(result[, c("emp_mi_github", "theory_mi_github", "theory_mi_ours_median_spikes", "theory_mi_ours_fixed_spikes", "theory_mi_ours_free_spikes")], na.rm = TRUE)
plot(result$umi, result$emp_mi_github, type = "b", log = "x", pch = 19, ylim = ylim,
     xlab = "Current-view UMI target", ylab = "MI", main = "GSE123025 cell-side overlap MI")
lines(result$umi, result$theory_mi_github, type = "b", pch = 19, col = "#D55E00")
lines(result$umi, result$theory_mi_ours_median_spikes, type = "b", pch = 19, col = "#CC79A7")
lines(result$umi, result$theory_mi_ours_fixed_spikes, type = "b", pch = 19, col = "#0072B2")
lines(result$umi, result$theory_mi_ours_free_spikes, type = "b", pch = 19, col = "#009E73")
legend("topleft",
       legend = c("Empirical overlap MI", "Theory, GitHub median MP/spikes", "Theory, ours median MP/spikes", "Theory, ours fixed-gamma quantile spikes", "Theory, ours free-gamma quantile spikes"),
       col = c("black", "#D55E00", "#CC79A7", "#0072B2", "#009E73"), lty = 1, pch = 19, bty = "n")
dev.off()

html <- c(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>MI comparison</title>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:28px;max-width:1120px} table{border-collapse:collapse;font-size:13px} th,td{border:1px solid #ddd;padding:6px 8px;text-align:right} th{background:#f4f6f8} img{max-width:100%;border:1px solid #ddd}</style>",
  "</head><body>",
  "<h1>GSE123025 overlap-MI comparison</h1>",
  sprintf("<p>Same disjoint molecular views: %d cells, %d HVGs, reference UMI target %s, top r=%d.</p>", n_cells, n_genes, u_ref, r_mi),
  "<h2>Summary</h2>",
  paste(capture.output(print(result, row.names = FALSE)), collapse = "<br>"),
  "<h2>Curves</h2><img src=\"mi_comparison_curves.png\">",
  "</body></html>"
)
writeLines(html, file.path(out_dir, "index.html"))

message("Wrote MI comparison to ", out_dir)
print(result)
