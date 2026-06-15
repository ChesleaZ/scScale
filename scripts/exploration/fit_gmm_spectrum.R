library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  eigenvalue_file <- args[[1]]
  out_dir <- args[[2]]
  dataset_label <- args[[3]]
  n_cells_use <- if (length(args) >= 4 && nzchar(args[[4]])) as.integer(args[[4]]) else NULL
} else {
  eigenvalue_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
  out_dir <- "outputs/exploration/gmm_spectrum"
  dataset_label <- "Jurkat"
  n_cells_use <- NULL
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_gaussian_modes <- 2
outlier_tail_quantile <- 0.98
max_iter <- 500
tol <- 1e-8
variance_floor <- 1e-10

fit_gmm_outlier_1d <- function(x, k, outlier_tail_quantile = 0.98, max_iter = 500, tol = 1e-8) {
  x <- x[is.finite(x)]
  n <- length(x)
  outlier_lower <- as.numeric(quantile(x, outlier_tail_quantile, names = FALSE))
  outlier_density <- 1 / (max(x) - outlier_lower)
  fit_x <- x[x < outlier_lower]
  is_outlier <- x >= outlier_lower
  is_bulk <- !is_outlier
  starts <- as.numeric(quantile(fit_x, seq(0.25, 0.75, length.out = k), names = FALSE))
  means <- starts
  variances <- rep(var(fit_x), k) / k
  outlier_weight <- mean(is_outlier)
  weights <- c(rep((1 - outlier_weight) / k, k), outlier_weight)
  resp <- cbind(
    sapply(seq_len(k), function(j) as.numeric(is_bulk) / k),
    as.numeric(is_outlier)
  )
  resp <- resp / rowSums(resp)
  loglik_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    gaussian_dens <- sapply(seq_len(k), function(j) {
      weights[j] * pmax(dnorm(x, mean = means[j], sd = sqrt(variances[j])), .Machine$double.xmin)
    })
    gaussian_dens[is_outlier, ] <- 0
    outlier_dens <- weights[k + 1] * ifelse(is_outlier, outlier_density, 0)
    dens <- cbind(gaussian_dens, outlier_dens)
    total <- pmax(rowSums(dens), .Machine$double.xmin)
    resp <- dens / total
    resp[is_outlier, seq_len(k)] <- 0
    resp[is_outlier, k + 1] <- 1
    resp[is_bulk, k + 1] <- 0
    resp[is_bulk, ] <- resp[is_bulk, ] / rowSums(resp[is_bulk, ])
    loglik_trace[iter] <- sum(log(total))

    nk <- colSums(resp[, seq_len(k), drop = FALSE])
    weights <- colSums(resp) / n
    means <- colSums(resp[, seq_len(k), drop = FALSE] * x) / nk
    variances <- colSums(resp[, seq_len(k), drop = FALSE] * (x - rep(means, each = n))^2) / nk
    variances <- pmax(variances, variance_floor)

    ord <- order(means)
    weights <- c(weights[ord], weights[k + 1])
    means <- means[ord]
    variances <- variances[ord]
    resp <- cbind(resp[, ord, drop = FALSE], resp[, k + 1])

    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      break
    }
  }

  list(
    parameters = data.frame(
      component = c(paste0("g", seq_len(k)), "outlier"),
      weight = weights,
      p_hat = weights * n,
      mean = c(means, NA),
      variance = c(variances, NA),
      sd = c(sqrt(variances), NA),
      outlier_lower = outlier_lower,
      n = n,
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    resp = resp,
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace)
  )
}

make_gmm_diagnostics <- function(x, fit, prefix, x_label, transform_inverse = identity) {
  params <- fit$parameters
  k <- nrow(params)
  x_grid <- seq(min(x), max(x), length.out = 20000)
  outlier_lower <- params$outlier_lower[1]
  outlier_density <- 1 / (max(x) - outlier_lower)

  curve_df <- do.call(rbind, lapply(seq_len(k), function(j) {
    density <- if (params$component[j] == "outlier") {
      ifelse(x_grid >= outlier_lower, params$weight[j] * outlier_density, 0)
    } else {
      params$weight[j] * dnorm(x_grid, mean = params$mean[j], sd = params$sd[j])
    }
    data.frame(
      component = params$component[j],
      x = x_grid,
      density = density
    )
  }))
  mix_df <- aggregate(density ~ x, curve_df, sum)
  dx <- c(diff(mix_df$x), tail(diff(mix_df$x), 1))
  mix_df$cdf <- cumsum(mix_df$density * dx)
  mix_df$cdf <- mix_df$cdf / max(mix_df$cdf)

  empirical_df <- data.frame(
    x = sort(x),
    empirical_cdf = seq_along(x) / length(x)
  )
  empirical_df$model_cdf <- approx(mix_df$x, mix_df$cdf, xout = empirical_df$x, rule = 2)$y

  qq_probs <- seq(0.0025, 0.9975, length.out = 600)
  qq_df <- data.frame(
    prob = qq_probs,
    empirical = as.numeric(quantile(x, qq_probs, names = FALSE)),
    model = approx(mix_df$cdf, mix_df$x, xout = qq_probs, ties = "ordered", rule = 2)$y
  )

  write.csv(params, file.path(out_dir, paste0(prefix, "_parameters.csv")), row.names = FALSE)
  write.csv(fit$loglik_trace, file.path(out_dir, paste0(prefix, "_loglik.csv")), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, paste0(prefix, "_empirical_vs_model_cdf.csv")), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, paste0(prefix, "_qq.csv")), row.names = FALSE)

  density_plot <- ggplot(data.frame(x = x), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = curve_df, aes(x, density, color = component), linewidth = 0.8) +
    geom_line(data = mix_df, aes(x, density), color = "black", linewidth = 1.1) +
    labs(
      title = paste0(prefix, " two-mode GMM + outlier density"),
      subtitle = paste0("Black curve is total density; outlier starts at ", signif(outlier_lower, 3)),
      x = x_label,
      y = "Density",
      color = "Component"
    )

  cdf_plot <- ggplot(empirical_df, aes(x)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    labs(
      title = paste0(prefix, " two-mode GMM + outlier CDF"),
      subtitle = "Grey: empirical CDF; blue: fitted model CDF",
      x = x_label,
      y = "CDF"
    )

  qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.35, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = qq_lim, ylim = qq_lim) +
    labs(
      title = paste0(prefix, " two-mode GMM + outlier Q-Q"),
      x = "Fitted model quantile",
      y = "Empirical quantile"
    )

  ggsave(file.path(out_dir, paste0(prefix, "_density.png")), density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, paste0(prefix, "_cdf.png")), cdf_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, paste0(prefix, "_qq.png")), qq_plot, width = 6, height = 6, dpi = 220)

  invisible(list(curve = curve_df, mixture = mix_df, empirical = empirical_df, qq = qq_df))
}

eig_df <- read.csv(eigenvalue_file)
n_cells <- if (is.null(n_cells_use)) max(eig_df$n_cells) else n_cells_use
eigvals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells], decreasing = TRUE)
eigvals <- eigvals[is.finite(eigvals) & eigvals > 0]

raw_fit <- fit_gmm_outlier_1d(eigvals, n_gaussian_modes, outlier_tail_quantile, max_iter = max_iter, tol = tol)
log_fit <- fit_gmm_outlier_1d(log10(eigvals), n_gaussian_modes, outlier_tail_quantile, max_iter = max_iter, tol = tol)

raw_assign <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  raw_component = raw_fit$parameters$component[max.col(raw_fit$resp, ties.method = "first")]
)
raw_assign <- cbind(raw_assign, setNames(as.data.frame(raw_fit$resp), paste0("raw_", raw_fit$parameters$component)))

log_assign <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  log10_eigenvalue = log10(eigvals),
  log_component = log_fit$parameters$component[max.col(log_fit$resp, ties.method = "first")]
)
log_assign <- cbind(log_assign, setNames(as.data.frame(log_fit$resp), paste0("log_", log_fit$parameters$component)))

write.csv(raw_assign, file.path(out_dir, "raw_gmm_assignments.csv"), row.names = FALSE)
write.csv(log_assign, file.path(out_dir, "log10_gmm_assignments.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

make_gmm_diagnostics(eigvals, raw_fit, "raw", "Covariance eigenvalue")
make_gmm_diagnostics(log10(eigvals), log_fit, "log10", "log10 covariance eigenvalue")

scree_plot <- ggplot(log_assign, aes(rank, eigenvalue, color = log_component)) +
  geom_point(size = 0.8, alpha = 0.85) +
  scale_y_log10() +
  labs(
    title = paste0(dataset_label, " two-mode GMM + outlier spectrum assignments"),
    subtitle = "MAP assignment from model fit on log10 eigenvalues",
    x = "Eigenvalue rank",
    y = "Covariance eigenvalue",
    color = "Component"
  )

ggsave(file.path(out_dir, "log10_gmm_scree_assignments.png"), scree_plot, width = 8, height = 5.5, dpi = 220)

print(raw_fit$parameters)
print(log_fit$parameters)
message("Dataset: ", dataset_label)
message("Eigenvalue file: ", eigenvalue_file)
message("n_cells: ", n_cells)
message("Wrote two-mode GMM plus outlier spectrum outputs to ", out_dir)
