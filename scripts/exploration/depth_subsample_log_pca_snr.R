#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: depth_subsample_log_pca_snr.R <counts_csv> <out_dir> [dataset_name]", call. = FALSE)
}

counts_csv <- args[[1]]
out_dir <- args[[2]]
dataset_name <- ifelse(length(args) >= 3, args[[3]], "dataset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

set.seed(as.integer(Sys.getenv("DEPTH_SEED", "1")))
rates <- as.numeric(strsplit(Sys.getenv("DEPTH_RATES", "0.05,0.1,0.2,0.5,1"), ",")[[1]])
fit_probs <- seq(0.05, 0.80, by = 0.025)
grid_n <- 20000

message("Reading counts from ", counts_csv)
counts_df <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
Y_full <- as.matrix(counts_df)
storage.mode(Y_full) <- "integer"

thin_counts <- function(Y, rate) {
  if (rate >= 1) {
    return(Y)
  }
  matrix(rbinom(length(Y), size = as.vector(Y), prob = rate), nrow = nrow(Y), ncol = ncol(Y))
}

umi_peff_summary <- function(Y) {
  umi <- colSums(Y)
  sum_sq <- colSums(Y^2)
  p_eff <- ifelse(sum_sq > 0, umi^2 / sum_sq, NA_real_)
  effective_depth <- ifelse(p_eff > 0, umi / p_eff, NA_real_)
  data.frame(
    mean_umi = mean(umi),
    median_umi = median(umi),
    mean_p_eff = mean(p_eff, na.rm = TRUE),
    median_p_eff = median(p_eff, na.rm = TRUE),
    mean_umi_per_p_eff = mean(effective_depth, na.rm = TRUE),
    median_umi_per_p_eff = median(effective_depth, na.rm = TRUE)
  )
}

log_normalized <- function(Y) {
  lib <- colSums(Y)
  lib[lib <= 0] <- 1
  log1p(sweep(Y, 2, lib, "/") * 1e4)
}

cov_eigenvalues <- function(X) {
  Xc <- sweep(X, 1, rowMeans(X), "-")
  eig <- svd(Xc, nu = 0, nv = 0)$d^2 / (ncol(Xc) - 1)
  sort(eig, decreasing = TRUE)
}

mp_density <- function(x, gamma, sigma2) {
  if (!is.finite(gamma) || !is.finite(sigma2) || gamma <= 0 || sigma2 <= 0) {
    return(rep(0, length(x)))
  }
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_positive_mass <- function(gamma) {
  ifelse(gamma > 1, 1 / gamma, 1)
}

mp_positive_density <- function(x, gamma, sigma2) {
  mp_density(x, gamma, sigma2) / mp_positive_mass(gamma)
}

mp_positive_quantiles <- function(gamma, sigma2, probs) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 0) + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_positive_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

fit_free_gamma <- function(vals, gamma_start) {
  vals <- vals[is.finite(vals) & vals > 1e-10]
  empirical_q <- as.numeric(quantile(vals, fit_probs, names = FALSE))
  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_positive_quantiles(gamma, sigma2, fit_probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }
  starts <- expand.grid(
    gamma = unique(pmax(c(gamma_start, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8), 1e-4)),
    sigma2 = unique(pmax(c(median(vals), mean(vals), quantile(vals, 0.75, names = FALSE)), .Machine$double.eps))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = log(c(1e-4, min(vals) / 100)),
      upper = log(c(50, max(vals) * 100)),
      control = list(maxit = 1500)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma <- unname(exp(best$par[1]))
  sigma2 <- unname(exp(best$par[2]))
  c(
    gamma = gamma,
    sigma2 = sigma2,
    mp_lower = sigma2 * (1 - sqrt(gamma))^2,
    mp_upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = best$value,
    convergence = best$convergence
  )
}

fit_bulk_iterative <- function(vals, gamma_theory, max_iter = 4) {
  positive <- vals[is.finite(vals) & vals > 1e-10]
  fit_vals <- positive
  fit <- NULL
  for (iter in seq_len(max_iter)) {
    fit <- fit_free_gamma(fit_vals, gamma_theory)
    next_vals <- positive[positive <= fit["mp_upper"]]
    if (length(next_vals) < 50 || length(next_vals) == length(fit_vals)) {
      break
    }
    fit_vals <- next_vals
  }
  c(fit, n_fit = length(fit_vals))
}

eig_rows <- list()
summary_rows <- list()
curve_rows <- list()

for (rate in rates) {
  message("Analyzing depth rate = ", rate)
  Y <- thin_counts(Y_full, rate)
  depth_summary <- umi_peff_summary(Y)
  eig <- cov_eigenvalues(log_normalized(Y))
  n_cells <- ncol(Y)
  n_genes <- nrow(Y)
  fit <- fit_bulk_iterative(eig, n_genes / n_cells)
  edge <- fit["mp_upper"]
  spikes <- eig[eig > edge]

  eig_rows[[length(eig_rows) + 1]] <- data.frame(
    dataset = dataset_name,
    depth_rate = rate,
    n_cells = n_cells,
    n_genes = n_genes,
    rank = seq_along(eig),
    eigenvalue = eig,
    log10_eigenvalue = log10(pmax(eig, .Machine$double.xmin))
  )

  summary_rows[[length(summary_rows) + 1]] <- cbind(
    data.frame(
      dataset = dataset_name,
      depth_rate = rate,
      n_cells = n_cells,
      n_genes = n_genes,
      gamma = fit["gamma"],
      sigma2 = fit["sigma2"],
      mp_lower = fit["mp_lower"],
      mp_upper = edge,
      objective = fit["objective"],
      convergence = fit["convergence"],
      n_fit = fit["n_fit"],
      lambda1 = eig[1],
      n_spikes = length(spikes),
      lambda1_over_edge = eig[1] / edge,
      spectral_snr = sqrt(eig[1] / edge),
      accumulated_snr = sqrt((n_cells - 1) * eig[1] / edge),
      normalized_spike_excess = sum(spikes / edge - 1)
    ),
    depth_summary
  )

  x <- seq(max(fit["mp_lower"], 1e-12), fit["mp_upper"], length.out = 1000)
  curve_rows[[length(curve_rows) + 1]] <- data.frame(
    depth_rate = rate,
    eigenvalue = x,
    density = mp_positive_density(x, fit["gamma"], fit["sigma2"])
  )
}

eig_df <- do.call(rbind, eig_rows)
summary_df <- do.call(rbind, summary_rows)
curve_df <- do.call(rbind, curve_rows)
rownames(summary_df) <- NULL

for (col in setdiff(names(summary_df), "dataset")) {
  summary_df[[col]] <- as.numeric(summary_df[[col]])
}

write.csv(eig_df, file.path(out_dir, paste0(dataset_name, "_depth_subsample_eigenvalues.csv")), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, paste0(dataset_name, "_depth_subsample_snr_summary.csv")), row.names = FALSE)
write.csv(curve_df, file.path(out_dir, paste0(dataset_name, "_depth_subsample_mp_curves.csv")), row.names = FALSE)

summary_df$depth_label <- factor(
  paste0("rate ", summary_df$depth_rate, "\nUMI/p_eff ", signif(summary_df$median_umi_per_p_eff, 3)),
  levels = paste0("rate ", summary_df$depth_rate, "\nUMI/p_eff ", signif(summary_df$median_umi_per_p_eff, 3))
)
eig_df$depth_label <- factor(
  paste0("rate ", eig_df$depth_rate),
  levels = paste0("rate ", summary_df$depth_rate)
)
curve_df$depth_label <- factor(
  paste0("rate ", curve_df$depth_rate),
  levels = paste0("rate ", summary_df$depth_rate)
)

theme_set(theme_classic(base_size = 13))

snr_long <- rbind(
  data.frame(depth_rate = summary_df$depth_rate, metric = "lambda1", value = summary_df$lambda1),
  data.frame(depth_rate = summary_df$depth_rate, metric = "MP upper edge", value = summary_df$mp_upper),
  data.frame(depth_rate = summary_df$depth_rate, metric = "spectral_snr", value = summary_df$spectral_snr),
  data.frame(depth_rate = summary_df$depth_rate, metric = "accumulated_snr", value = summary_df$accumulated_snr),
  data.frame(depth_rate = summary_df$depth_rate, metric = "median UMI / p_eff", value = summary_df$median_umi_per_p_eff)
)

snr_plot <- ggplot(snr_long, aes(depth_rate, value)) +
  geom_line(linewidth = 0.9, color = "#1F78B4") +
  geom_point(size = 2.1, color = "#1F78B4") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": depth subsampling SNR diagnostics"),
    x = "Sampling ratio",
    y = "Value (log scale)"
  )

bulk_plot_df <- do.call(rbind, lapply(seq_len(nrow(summary_df)), function(i) {
  rate <- summary_df$depth_rate[i]
  vals <- eig_df$eigenvalue[eig_df$depth_rate == rate]
  upper <- summary_df$mp_upper[i]
  data.frame(depth_rate = rate, eigenvalue = vals[vals <= upper * 1.35])
}))
bulk_plot_df$depth_label <- factor(
  paste0("rate ", bulk_plot_df$depth_rate),
  levels = paste0("rate ", summary_df$depth_rate)
)

bulk_plot <- ggplot(bulk_plot_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_line(data = curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.9) +
  geom_vline(data = summary_df, aes(xintercept = mp_upper), color = "#C43C39", linetype = "dashed") +
  facet_wrap(~ depth_label, scales = "free", ncol = 2) +
  labs(
    title = paste0(dataset_name, ": MP bulk fits by depth"),
    x = "Log-PCA covariance eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, paste0(dataset_name, "_depth_subsample_snr_diagnostics.png")), snr_plot, width = 10, height = 8, dpi = 220)
ggsave(file.path(out_dir, paste0(dataset_name, "_depth_subsample_mp_bulk_fits.png")), bulk_plot, width = 10, height = 8, dpi = 220)

print(summary_df)
message("Wrote depth subsampling outputs to ", out_dir)
