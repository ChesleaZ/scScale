library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/jurkat_exact_mp_pipeline"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_counts_use <- c(2000, 5000, 10000)
knee_max_rank <- 100

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

find_scree_knee <- function(eigvals, max_rank = 100) {
  eigvals <- sort(eigvals, decreasing = TRUE)
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

estimate_exact_mp_after_spikes <- function(eigvals, n, k_spikes) {
  eigvals <- sort(eigvals, decreasing = TRUE)
  p <- length(eigvals)
  noise_vals <- eigvals[(k_spikes + 1):p]
  noise_vals <- noise_vals[noise_vals > 1e-10]
  p_eff <- length(noise_vals)

  gamma <- p_eff / n
  sigma2 <- mean(noise_vals)
  mp_lower <- sigma2 * (1 - sqrt(gamma))^2
  mp_upper <- sigma2 * (1 + sqrt(gamma))^2

  data.frame(
    n_cells = n,
    n_genes = p,
    k_spikes = k_spikes,
    n_noise_positive = p_eff,
    gamma = gamma,
    sigma2 = sigma2,
    mp_lower = mp_lower,
    mp_upper = mp_upper,
    knee_max_rank = knee_max_rank
  )
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
summary_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells_summary.csv"))

cell_counts <- intersect(cell_counts_use, sort(unique(eig_df$n_cells)))

fit_df <- do.call(rbind, lapply(cell_counts, function(n_cells) {
  eigvals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  k_spikes <- find_scree_knee(eigvals, max_rank = knee_max_rank)
  fit <- estimate_exact_mp_after_spikes(eigvals, n = n_cells, k_spikes = k_spikes)
  row <- summary_df[summary_df$n_cells == n_cells, ]
  fit$trace_mean_all <- row$mean_eigenvalue
  fit$max_eigenvalue <- row$max_eigenvalue
  fit
}))

labels <- paste0("n = ", format(cell_counts, big.mark = ","))
names(labels) <- cell_counts

residual_df <- do.call(rbind, lapply(seq_len(nrow(fit_df)), function(i) {
  fit <- fit_df[i, ]
  vals <- sort(eig_df$eigenvalue[eig_df$n_cells == fit$n_cells], decreasing = TRUE)
  data.frame(
    n_cells = fit$n_cells,
    rank = seq_along(vals),
    eigenvalue = vals,
    component = ifelse(seq_along(vals) <= fit$k_spikes, "detected spike", "residual bulk")
  )
}))

curve_df <- do.call(rbind, lapply(seq_len(nrow(fit_df)), function(i) {
  fit <- fit_df[i, ]
  x <- seq(max(fit$mp_lower, 1e-12), fit$mp_upper, length.out = 1000)
  data.frame(
    n_cells = fit$n_cells,
    eigenvalue = x,
    density = mp_density(x, fit$gamma, fit$sigma2)
  )
}))

fit_df$n_cells_label <- factor(labels[as.character(fit_df$n_cells)], levels = labels)
residual_df$n_cells_label <- factor(labels[as.character(residual_df$n_cells)], levels = labels)
curve_df$n_cells_label <- factor(labels[as.character(curve_df$n_cells)], levels = labels)

write.csv(fit_df, file.path(out_dir, "jurkat_exact_mp_parameters.csv"), row.names = FALSE)
write.csv(residual_df, file.path(out_dir, "jurkat_exact_mp_classified_eigenvalues.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

bulk_plot <- ggplot(subset(residual_df, component == "residual bulk" & eigenvalue > 1e-10), aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 120,
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
  geom_vline(data = fit_df, aes(xintercept = mp_upper), color = "#C43C39", linetype = "dashed") +
  geom_text(
    data = fit_df,
    aes(
      x = Inf,
      y = Inf,
      label = paste0(
        "spikes: ", k_spikes,
        "\nsigma2: ", signif(sigma2, 3),
        "\ngamma: ", signif(gamma, 3)
      )
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.25,
    size = 3.1
  ) +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Jurkat HVG residual bulk vs exact MP",
    subtitle = "Spikes selected by scree elbow; MP scale is residual trace, not a shape fit",
    x = "Covariance eigenvalue",
    y = "Density"
  )

scree_plot <- ggplot(subset(residual_df, rank <= 80), aes(rank, eigenvalue, color = component)) +
  geom_line(color = "grey45") +
  geom_point(size = 1.8) +
  geom_hline(data = fit_df, aes(yintercept = mp_upper), color = "#C43C39", linetype = "dashed") +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 1) +
  scale_y_log10() +
  labs(
    title = "Jurkat exact-MP spike calls",
    x = "Eigenvalue rank",
    y = "Eigenvalue",
    color = NULL
  )

ggsave(file.path(out_dir, "jurkat_exact_mp_bulk_by_cells.png"), bulk_plot, width = 8, height = 10, dpi = 220)
ggsave(file.path(out_dir, "jurkat_exact_mp_scree_by_cells.png"), scree_plot, width = 8, height = 10, dpi = 220)

print(fit_df)
message("Wrote Jurkat exact MP pipeline outputs to ", out_dir)
