library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

set.seed(1)

data_dir <- "data/Jurkat/sample_filtered_feature_bc_matrix"
out_dir <- "outputs/exploration/jurkat_eigenvalues"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_variable_genes <- 2000
mp_bulk_quantile <- 0.90

mp_unit_quantiles <- function(gamma, probs, grid_n = 5000) {
  gamma <- min(max(gamma, 1e-4), 0.999)
  lower <- (1 - sqrt(gamma))^2
  upper <- (1 + sqrt(gamma))^2
  x <- seq(lower + 1e-12, upper - 1e-12, length.out = grid_n)
  density <- sqrt((upper - x) * (x - lower)) / (2 * pi * gamma * x)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(density * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

fit_mp_bulk <- function(eigvals, bulk_quantile = 0.90) {
  fit_vals <- eigvals[eigvals > 0 & eigvals <= quantile(eigvals, bulk_quantile)]
  probs <- seq(0.05, 0.95, by = 0.05)
  empirical_q <- as.numeric(quantile(fit_vals, probs))

  objective <- function(par) {
    gamma <- plogis(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- sigma2 * mp_unit_quantiles(gamma, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  start_gamma <- 0.85
  start_sigma2 <- median(fit_vals)
  fit <- optim(
    par = c(qlogis(start_gamma), log(start_sigma2)),
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = 1000)
  )

  gamma <- plogis(fit$par[1])
  sigma2 <- exp(fit$par[2])
  data.frame(
    gamma = gamma,
    sigma2 = sigma2,
    lower = sigma2 * (1 - sqrt(gamma))^2,
    upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = fit$value,
    n_fit = length(fit_vals),
    bulk_quantile = bulk_quantile
  )
}

message("Loading 10x matrix from ", data_dir)
counts <- Read10X(data.dir = data_dir)
seu <- CreateSeuratObject(counts)

message("Normalizing and selecting variable genes")
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = n_variable_genes,
  verbose = FALSE
)

genes_use <- VariableFeatures(seu)
X <- GetAssayData(seu, assay = DefaultAssay(seu), layer = "data")[genes_use, ]

p <- nrow(X)
n <- ncol(X)
message("Genes used: ", p)
message("Cells used: ", n)
message("Aspect ratio p/n: ", signif(p / n, 4))

message("Computing centered gene covariance")
gene_means <- Matrix::rowMeans(X)

# cov(X_centered) = (X X' - n * mean mean') / (n - 1)
cov_mat <- as.matrix(Matrix::tcrossprod(X))
cov_mat <- (cov_mat - n * tcrossprod(gene_means)) / (n - 1)
cov_mat <- (cov_mat + t(cov_mat)) / 2

message("Computing all covariance eigenvalues")
eigvals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
eigvals[eigvals < 0 & eigvals > -1e-8] <- 0
eigvals <- sort(eigvals, decreasing = TRUE)

message("Fitting empirical MP bulk")
mp_fit <- fit_mp_bulk(eigvals, bulk_quantile = mp_bulk_quantile)
mp_lower <- mp_fit$lower
mp_upper <- mp_fit$upper
n_above_mp <- sum(eigvals > mp_upper)

mp_curve <- data.frame(
  eigenvalue = seq(mp_lower, mp_upper, length.out = 1000)
)
mp_curve$density <- mp_density(mp_curve$eigenvalue, mp_fit$gamma, mp_fit$sigma2)
mp_curve$log10_eigenvalue <- log10(mp_curve$eigenvalue)
mp_curve$log10_density <- mp_curve$density * mp_curve$eigenvalue * log(10)

summary_df <- data.frame(
  metric = c(
    "n_genes",
    "n_cells",
    "p_over_n",
    "min_eigenvalue",
    "median_eigenvalue",
    "mean_eigenvalue",
    "max_eigenvalue",
    "mp_fit_bulk_quantile",
    "mp_fit_n_eigenvalues",
    "mp_fit_gamma",
    "mp_fit_sigma2",
    "mp_fit_lower",
    "mp_fit_upper",
    "mp_fit_objective",
    "n_above_mp"
  ),
  value = c(
    p,
    n,
    p / n,
    min(eigvals),
    median(eigvals),
    mean(eigvals),
    max(eigvals),
    mp_fit$bulk_quantile,
    mp_fit$n_fit,
    mp_fit$gamma,
    mp_fit$sigma2,
    mp_lower,
    mp_upper,
    mp_fit$objective,
    n_above_mp
  )
)

eig_df <- data.frame(
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  log10_eigenvalue = log10(pmax(eigvals, .Machine$double.eps)),
  above_mp = eigvals > mp_upper
)

write.csv(eig_df, file.path(out_dir, "jurkat_eigenvalues.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "jurkat_eigenvalue_summary.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 14))

scree <- ggplot(eig_df, aes(rank, eigenvalue)) +
  geom_line(linewidth = 0.5, color = "grey35") +
  geom_point(aes(color = above_mp), size = 0.7, alpha = 0.75) +
  geom_hline(yintercept = mp_upper, color = "#C43C39", linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "grey35", "TRUE" = "#1F78B4"), guide = "none") +
  labs(
    title = "Jurkat covariance spectrum",
    subtitle = paste0(n_above_mp, " eigenvalues above estimated MP upper edge"),
    x = "Eigenvalue rank",
    y = "Eigenvalue"
  )

log_scree <- ggplot(eig_df, aes(rank, eigenvalue)) +
  geom_line(linewidth = 0.5, color = "grey35") +
  geom_point(aes(color = above_mp), size = 0.7, alpha = 0.75) +
  geom_hline(yintercept = mp_upper, color = "#C43C39", linetype = "dashed") +
  scale_y_log10() +
  scale_color_manual(values = c("FALSE" = "grey35", "TRUE" = "#1F78B4"), guide = "none") +
  labs(
    title = "Jurkat covariance spectrum, log scale",
    x = "Eigenvalue rank",
    y = "Eigenvalue"
  )

hist_linear <- ggplot(eig_df, aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 80,
    fill = "grey70",
    color = "white",
    linewidth = 0.2
  ) +
  geom_line(data = mp_curve, aes(eigenvalue, density), color = "#1F78B4", linewidth = 1) +
  geom_vline(xintercept = mp_upper, color = "#C43C39", linetype = "dashed") +
  labs(
    title = "Eigenvalue distribution with fitted MP bulk",
    x = "Eigenvalue",
    y = "Density"
  )

hist_log <- ggplot(eig_df, aes(log10_eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 80,
    fill = "grey70",
    color = "white",
    linewidth = 0.2
  ) +
  geom_line(
    data = mp_curve,
    aes(log10_eigenvalue, log10_density),
    color = "#1F78B4",
    linewidth = 1
  ) +
  geom_vline(xintercept = log10(mp_upper), color = "#C43C39", linetype = "dashed") +
  labs(
    title = "Eigenvalue distribution with fitted MP bulk, log10 scale",
    x = "log10 eigenvalue",
    y = "Density"
  )

combined <- (scree | log_scree) / (hist_linear | hist_log)

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_distributions.png"),
  combined,
  width = 13,
  height = 9,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_scree.png"),
  scree,
  width = 8,
  height = 5,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_histogram_log10.png"),
  hist_log,
  width = 8,
  height = 5,
  dpi = 220
)

print(summary_df)
message("Wrote outputs to ", out_dir)
