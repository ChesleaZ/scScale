library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/two_mp_mixture_em"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cells_use <- NULL
max_iter <- 500
tol <- 1e-7
density_floor_multiplier <- 1e-10
spike_tail_quantile <- 0.95

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

weighted_fit_sigma2 <- function(vals, weights, gamma, density_floor) {
  positive_vals <- vals[vals > 0]
  lower <- log(max(min(positive_vals) * 0.05, .Machine$double.eps))
  upper <- log(max(positive_vals) * 20)

  objective <- function(log_sigma2) {
    sigma2 <- exp(log_sigma2)
    d <- pmax(mp_density(vals, gamma, sigma2), density_floor)
    -sum(weights * log(d))
  }

  opt <- optimize(objective, interval = c(lower, upper))
  exp(opt$minimum)
}

order_mp_components <- function(resp, counts, sigma2) {
  ord <- order(sigma2)
  list(
    resp = cbind(resp[, ord, drop = FALSE], resp[, 3]),
    counts = c(counts[ord], counts[3]),
    sigma2 = sigma2[ord]
  )
}

fit_two_mp_em <- function(vals, n_cells) {
  vals <- vals[is.finite(vals) & vals > 0]
  p <- length(vals)
  spike_lower <- as.numeric(quantile(vals, spike_tail_quantile))
  spike_density_value <- 1 / (max(vals) - spike_lower)
  density_floor <- density_floor_multiplier * spike_density_value

  p3_init <- max(20, sum(vals >= spike_lower))
  p2_init <- round(0.35 * (p - p3_init))
  p1_init <- p - p2_init - p3_init
  counts <- c(p1_init, p2_init, p3_init)

  sigma2 <- c(
    median(vals[vals <= quantile(vals, 0.55)]),
    median(vals[vals > quantile(vals, 0.55) & vals <= quantile(vals, 0.95)])
  )

  resp <- matrix(1 / 3, nrow = p, ncol = 3)
  loglik_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    gamma <- counts[1:2] / n_cells
    pi <- counts / p

    d1 <- pmax(mp_density(vals, gamma[1], sigma2[1]), density_floor)
    d2 <- pmax(mp_density(vals, gamma[2], sigma2[2]), density_floor)
    d3 <- ifelse(vals >= spike_lower, spike_density_value, density_floor)

    weighted_density <- cbind(pi[1] * d1, pi[2] * d2, pi[3] * d3)
    total_density <- rowSums(weighted_density)
    resp <- weighted_density / total_density

    loglik <- sum(log(total_density))
    loglik_trace[iter] <- loglik

    counts <- pmax(colSums(resp), c(5, 5, 1))
    counts <- counts * p / sum(counts)
    gamma <- counts[1:2] / n_cells

    sigma2 <- c(
      weighted_fit_sigma2(vals, resp[, 1], gamma[1], density_floor),
      weighted_fit_sigma2(vals, resp[, 2], gamma[2], density_floor)
    )

    ordered <- order_mp_components(resp, counts, sigma2)
    resp <- ordered$resp
    counts <- ordered$counts
    sigma2 <- ordered$sigma2

    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      break
    }
  }

  gamma <- counts[1:2] / n_cells
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2

  list(
    resp = resp,
    parameters = data.frame(
      component = c("noise_mp", "intrinsic_mp", "spike_outlier"),
      p_hat = counts,
      weight = counts / p,
      gamma = c(gamma, NA),
      sigma2 = c(sigma2, NA),
      mp_lower = c(lower, NA),
      mp_upper = c(upper, NA),
      n_cells = n_cells,
      n_eigenvalues = p,
      spike_lower = spike_lower,
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace)
  )
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
n_cells <- if (is.null(n_cells_use)) max(eig_df$n_cells) else n_cells_use
eigvals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells], decreasing = TRUE)
fit <- fit_two_mp_em(eigvals, n_cells)

assign_df <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  noise_mp = fit$resp[, 1],
  intrinsic_mp = fit$resp[, 2],
  spike_outlier = fit$resp[, 3]
)
assign_df$component <- c("noise_mp", "intrinsic_mp", "spike_outlier")[max.col(fit$resp, ties.method = "first")]

write.csv(fit$parameters, file.path(out_dir, "two_mp_mixture_em_parameters.csv"), row.names = FALSE)
write.csv(assign_df, file.path(out_dir, "two_mp_mixture_em_assignments.csv"), row.names = FALSE)
write.csv(fit$loglik_trace, file.path(out_dir, "two_mp_mixture_em_loglik.csv"), row.names = FALSE)

x_max <- as.numeric(quantile(eigvals, 0.995))
x <- seq(min(eigvals), x_max, length.out = 2000)
params <- fit$parameters
spike_lower <- params$spike_lower[1]
spike_density_value <- 1 / (max(eigvals) - spike_lower)

curve_df <- rbind(
  data.frame(
    component = "noise_mp",
    eigenvalue = x,
    density = params$weight[1] * mp_density(x, params$gamma[1], params$sigma2[1])
  ),
  data.frame(
    component = "intrinsic_mp",
    eigenvalue = x,
    density = params$weight[2] * mp_density(x, params$gamma[2], params$sigma2[2])
  ),
  data.frame(
    component = "spike_outlier",
    eigenvalue = x,
    density = ifelse(x >= spike_lower, params$weight[3] * spike_density_value, 0)
  )
)
mix_df <- aggregate(density ~ eigenvalue, curve_df, sum)

theme_set(theme_classic(base_size = 13))

density_plot <- ggplot(data.frame(eigenvalue = eigvals), aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 140, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_line(data = curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
  geom_line(data = mix_df, aes(eigenvalue, density), linewidth = 1.1, color = "black") +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(
    title = "Two-MP mixture EM fit",
    subtitle = paste0("Black curve is total mixture; spike density is restricted to eigenvalues >= ", signif(spike_lower, 3)),
    x = "Covariance eigenvalue",
    y = "Density",
    color = "Component"
  )

scree_plot <- ggplot(assign_df, aes(rank, eigenvalue, color = component)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Two-MP mixture EM assignments",
    x = "Eigenvalue rank",
    y = "Covariance eigenvalue",
    color = "MAP component"
  )

loglik_plot <- ggplot(fit$loglik_trace, aes(iter, loglik)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title = "EM log likelihood",
    x = "Iteration",
    y = "Log likelihood"
  )

ggsave(file.path(out_dir, "two_mp_mixture_em_density_fit.png"), density_plot, width = 9, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_mixture_em_scree_assignments.png"), scree_plot, width = 8, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_mixture_em_loglik.png"), loglik_plot, width = 7, height = 4.5, dpi = 220)

print(fit$parameters)
message("Wrote two-MP mixture EM outputs to ", out_dir)
