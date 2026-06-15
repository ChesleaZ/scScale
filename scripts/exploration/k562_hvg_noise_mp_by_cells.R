library(Seurat)
library(Matrix)
library(ggplot2)

set.seed(1)

dataset_name <- "k562_320k"
data_dir <- "data/K562/sample_filtered_feature_bc_matrix"
out_dir <- file.path("outputs", "k562_320k_hvg_noise_mp_by_cells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_variable_genes <- 2000
cell_counts_requested <- c(500, 1000, 2000, 5000, 10000, 25000, 50000, 100000, 200000, Inf)
knee_max_rank <- 100
n_bins <- 160

eig_file <- file.path(out_dir, "k562_hvg_eigenvalues_by_cells.csv")
summary_file <- file.path(out_dir, "k562_hvg_eigenvalues_by_cells_summary.csv")
noise_file <- file.path(out_dir, "k562_hvg_noise_eigenvalues_by_cells.csv")
fit_file <- file.path(out_dir, "k562_hvg_noise_threshold_mp_fit_by_cells.csv")
progress_file <- file.path(out_dir, "progress.log")

log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = progress_file, append = TRUE)
}

compute_cov_eigenvalues <- function(X) {
  n <- ncol(X)
  gene_means <- Matrix::rowMeans(X)
  cov_mat <- as.matrix(Matrix::tcrossprod(X))
  cov_mat <- (cov_mat - n * tcrossprod(gene_means)) / (n - 1)
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  eigvals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eigvals[eigvals < 0 & eigvals > -1e-8] <- 0
  sort(eigvals, decreasing = TRUE)
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

mp_quantiles <- function(gamma, sigma2, probs, grid_n = 6000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(lower + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  approx(cdf, x, xout = probs, ties = "ordered", rule = 2)$y
}

find_scree_knee <- function(eigvals, max_rank = 100) {
  k <- min(max_rank, length(eigvals) - 1)
  y <- log10(pmax(eigvals[seq_len(k)], .Machine$double.eps))
  x <- seq_along(y)
  x <- (x - min(x)) / (max(x) - min(x))
  y <- (y - min(y)) / (max(y) - min(y))

  x1 <- x[1]
  y1 <- y[1]
  x2 <- tail(x, 1)
  y2 <- tail(y, 1)

  dist <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  which.max(dist)
}

fit_mp <- function(vals) {
  fit_vals <- vals[vals > 1e-10]
  probs <- seq(0.02, 0.98, by = 0.02)
  empirical_q <- as.numeric(quantile(fit_vals, probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_quantiles(gamma, sigma2, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = c(0.2, 0.5, 1, 2, 5, 10),
    sigma2 = c(median(fit_vals), mean(fit_vals), quantile(fit_vals, 0.75))
  )

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$gamma[i], starts$sigma2[i])),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1500)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  gamma <- exp(best$par[1])
  sigma2 <- exp(best$par[2])

  data.frame(
    mp_gamma = gamma,
    mp_sigma2 = sigma2,
    mp_lower = sigma2 * (1 - sqrt(gamma))^2,
    mp_upper = sigma2 * (1 + sqrt(gamma))^2,
    mp_objective = best$value,
    n_positive_fit = length(fit_vals)
  )
}

if (file.exists(eig_file)) {
  log_step("Loading existing eigenvalue checkpoint from ", eig_file)
  eig_df <- read.csv(eig_file)
} else {
  eig_df <- data.frame()
}

if (!file.exists(eig_file) || length(setdiff(cell_counts_requested, unique(eig_df$n_cells))) > 0) {
  log_step("Loading 10x matrix from ", data_dir)
  counts <- Read10X(data.dir = data_dir)
  if (is.list(counts)) {
    log_step("Read10X returned feature types: ", paste(names(counts), collapse = ", "))
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
      log_step("Gene Expression feature type not found; using first matrix")
    }
  }
  seu <- CreateSeuratObject(counts)
  rm(counts)
  gc()

  log_step("Normalizing and selecting ", n_variable_genes, " HVGs")
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
  X_full <- GetAssayData(seu, assay = DefaultAssay(seu), layer = "data")[genes_use, ]
  n_cells_total <- ncol(X_full)
  cell_counts <- pmin(cell_counts_requested, n_cells_total)
  cell_counts <- unique(as.integer(cell_counts))

  log_step("Genes used: ", nrow(X_full))
  log_step("Total cells available: ", n_cells_total)
  log_step("Cell counts analyzed: ", paste(cell_counts, collapse = ", "))

  done_counts <- unique(eig_df$n_cells)
  for (n_cells in cell_counts) {
    if (n_cells %in% done_counts) {
      log_step("Skipping completed n_cells = ", n_cells)
      next
    }

    log_step("Computing eigenvalues for n_cells = ", n_cells)
    cell_idx <- if (n_cells == n_cells_total) {
      seq_len(n_cells_total)
    } else {
      sort(sample.int(n_cells_total, n_cells))
    }

    eigvals <- compute_cov_eigenvalues(X_full[, cell_idx, drop = FALSE])
    eig_new <- data.frame(
      n_cells = n_cells,
      rank = seq_along(eigvals),
      eigenvalue = eigvals,
      log10_eigenvalue = log10(pmax(eigvals, .Machine$double.eps))
    )

    eig_df <- rbind(eig_df, eig_new)
    write.csv(eig_df, eig_file, row.names = FALSE)
    log_step("Checkpointed eigenvalues through n_cells = ", n_cells)
    gc()
  }
} else {
  log_step("All requested spectra already present; skipping matrix load")
}

cell_counts <- sort(unique(eig_df$n_cells))
summary_df <- do.call(rbind, lapply(cell_counts, function(n_cells) {
  eigvals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  data.frame(
    n_cells = n_cells,
    n_genes = length(eigvals),
    max_eigenvalue = max(eigvals),
    median_eigenvalue = median(eigvals),
    mean_eigenvalue = mean(eigvals),
    n_eigenvalues_gt_0.5 = sum(eigvals > 0.5),
    n_eigenvalues_gt_2 = sum(eigvals > 2)
  )
}))
write.csv(summary_df, summary_file, row.names = FALSE)

noise_list <- vector("list", length(cell_counts))
fit_list <- vector("list", length(cell_counts))
curve_list <- vector("list", length(cell_counts))

for (i in seq_along(cell_counts)) {
  n_cells <- cell_counts[i]
  eigvals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  eigvals <- sort(eigvals, decreasing = TRUE)

  spike_count <- find_scree_knee(eigvals, max_rank = knee_max_rank)
  noise_upper <- eigvals[spike_count + 1]
  noise_vals <- eigvals[(spike_count + 1):length(eigvals)]
  fit <- fit_mp(noise_vals)

  x_max <- min(max(noise_vals), max(fit$mp_upper, noise_upper))
  x <- seq(max(fit$mp_lower, 1e-12), x_max, length.out = 1000)
  positive_mass <- mean(noise_vals > 1e-10)

  noise_list[[i]] <- data.frame(
    n_cells = n_cells,
    rank = seq_along(noise_vals) + spike_count,
    eigenvalue = noise_vals
  )
  fit_list[[i]] <- cbind(
    data.frame(
      n_cells = n_cells,
      spike_count = spike_count,
      noise_upper = noise_upper,
      n_noise = length(noise_vals),
      n_noise_zero = sum(noise_vals <= 1e-10)
    ),
    fit
  )
  curve_list[[i]] <- data.frame(
    n_cells = n_cells,
    eigenvalue = x,
    density = mp_density(x, fit$mp_gamma, fit$mp_sigma2) * positive_mass
  )
}

noise_df <- do.call(rbind, noise_list)
fit_df <- do.call(rbind, fit_list)
curve_df <- do.call(rbind, curve_list)

labels <- paste0("n = ", format(cell_counts, big.mark = ","))
names(labels) <- cell_counts
noise_df$n_cells_label <- factor(labels[as.character(noise_df$n_cells)], levels = labels)
fit_df$n_cells_label <- factor(labels[as.character(fit_df$n_cells)], levels = labels)
curve_df$n_cells_label <- factor(labels[as.character(curve_df$n_cells)], levels = labels)

write.csv(noise_df, noise_file, row.names = FALSE)
write.csv(fit_df, fit_file, row.names = FALSE)

theme_set(theme_classic(base_size = 13))

hist_plot <- ggplot(noise_df, aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = n_bins,
    fill = "grey72",
    color = "white",
    linewidth = 0.12
  ) +
  geom_line(
    data = curve_df,
    aes(eigenvalue, density),
    color = "#1F78B4",
    linewidth = 0.9
  ) +
  geom_vline(
    data = fit_df,
    aes(xintercept = noise_upper),
    color = "#C43C39",
    linetype = "dashed"
  ) +
  geom_text(
    data = fit_df,
    aes(
      x = Inf,
      y = Inf,
      label = paste0("spikes: ", spike_count, "\nthreshold: ", signif(noise_upper, 3))
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.25,
    size = 3.1
  ) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  coord_cartesian(xlim = c(0, max(fit_df$noise_upper) * 1.05)) +
  labs(
    title = "K562 HVG noise eigenvalue distributions with free MP overlays",
    subtitle = paste0("2,000 HVGs; noise threshold from scree elbow among top ", knee_max_rank, " ranks"),
    x = "Covariance eigenvalue",
    y = "Density"
  )

threshold_plot <- ggplot(fit_df, aes(n_cells, spike_count)) +
  geom_line(linewidth = 0.9, color = "grey25") +
  geom_point(size = 2.2, color = "#1F78B4") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "K562 automatically selected spike count",
    x = "Number of cells",
    y = "Spike count"
  )

ggsave(
  file.path(out_dir, "k562_hvg_noise_distributions_by_cells_free_mp.png"),
  hist_plot,
  width = 11,
  height = 14,
  dpi = 220
)

ggsave(
  file.path(out_dir, "k562_hvg_noise_spike_count_by_cells.png"),
  threshold_plot,
  width = 8,
  height = 5,
  dpi = 220
)

print(summary_df)
print(fit_df)
log_step("Done")
