library(ggplot2)

set.seed(7)

out_dir <- "outputs/exploration/simulated_spiked_mp_pipeline"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- 1000
n <- 5000
sigma2_true <- 0.35
spike_strengths <- c(4.0, 2.5, 1.5, 0.9, 0.55)

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

estimate_spikes_exact_mp <- function(eigvals, n, edge_buffer = 1.01, max_iter = 100) {
  eigvals <- sort(eigvals, decreasing = TRUE)
  p <- length(eigvals)
  k <- 0

  for (iter in seq_len(max_iter)) {
    noise_vals <- eigvals[(k + 1):p]
    p_eff <- length(noise_vals)
    gamma <- p_eff / n
    sigma2 <- mean(noise_vals)
    mp_lower <- sigma2 * (1 - sqrt(gamma))^2
    mp_upper <- sigma2 * (1 + sqrt(gamma))^2
    k_new <- sum(eigvals > edge_buffer * mp_upper)

    if (k_new == k) {
      return(data.frame(
        k = k,
        sigma2 = sigma2,
        gamma = gamma,
        mp_lower = mp_lower,
        mp_upper = mp_upper,
        edge_buffer = edge_buffer,
        n_noise = p_eff,
        iterations = iter
      ))
    }

    k <- min(k_new, p - 1)
  }

  stop("Spike-count iteration did not converge.")
}

simulate_spiked_covariance <- function(p, n, sigma2, spike_strengths) {
  r <- length(spike_strengths)
  raw_u <- matrix(rnorm(p * r), p, r)
  u <- qr.Q(qr(raw_u))

  factors <- matrix(rnorm(r * n), r, n)
  signal <- u %*% (sqrt(spike_strengths) * factors)
  noise <- matrix(rnorm(p * n, sd = sqrt(sigma2)), p, n)
  x <- signal + noise
  x <- x - rowMeans(x)

  cov_mat <- tcrossprod(x) / (n - 1)
  eigvals <- sort(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)

  list(eigvals = eigvals)
}

sim <- simulate_spiked_covariance(p, n, sigma2_true, spike_strengths)
fit <- estimate_spikes_exact_mp(sim$eigvals, n = n)

truth <- data.frame(
  p = p,
  n = n,
  gamma_true = p / n,
  sigma2_true = sigma2_true,
  n_spikes_true = length(spike_strengths),
  population_spike_eigenvalues = paste(signif(sigma2_true + spike_strengths, 4), collapse = ", ")
)

summary_df <- cbind(truth, fit)
write.csv(summary_df, file.path(out_dir, "spiked_mp_recovery_summary.csv"), row.names = FALSE)

eig_df <- data.frame(
  rank = seq_along(sim$eigvals),
  eigenvalue = sim$eigvals,
  component = ifelse(seq_along(sim$eigvals) <= fit$k, "detected spike", "residual bulk")
)
write.csv(eig_df, file.path(out_dir, "spiked_mp_eigenvalues.csv"), row.names = FALSE)

curve_x <- seq(fit$mp_lower, fit$mp_upper, length.out = 1000)
curve_df <- data.frame(
  eigenvalue = curve_x,
  density = mp_density(curve_x, fit$gamma, fit$sigma2)
)

theme_set(theme_classic(base_size = 13))

hist_plot <- ggplot(subset(eig_df, component == "residual bulk"), aes(eigenvalue)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 120,
    fill = "grey72",
    color = "white",
    linewidth = 0.12
  ) +
  geom_line(data = curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 1) +
  geom_vline(xintercept = fit$mp_upper, color = "#C43C39", linetype = "dashed") +
  labs(
    title = "Simulated spiked covariance: residual bulk vs exact MP",
    subtitle = paste0(
      "true spikes = ", length(spike_strengths),
      "; detected spikes = ", fit$k,
      "; true sigma2 = ", sigma2_true,
      "; estimated sigma2 = ", signif(fit$sigma2, 4)
    ),
    x = "Covariance eigenvalue",
    y = "Density"
  )

scree_plot <- ggplot(eig_df[seq_len(60), ], aes(rank, eigenvalue, color = component)) +
  geom_point(size = 2) +
  geom_line(color = "grey45") +
  geom_hline(yintercept = fit$mp_upper, color = "#C43C39", linetype = "dashed") +
  scale_y_log10() +
  labs(
    title = "Detected spikes in simulated covariance spectrum",
    x = "Eigenvalue rank",
    y = "Eigenvalue",
    color = NULL
  )

ggsave(file.path(out_dir, "simulated_spiked_mp_bulk.png"), hist_plot, width = 8, height = 5, dpi = 220)
ggsave(file.path(out_dir, "simulated_spiked_mp_scree.png"), scree_plot, width = 8, height = 5, dpi = 220)

print(summary_df)
message("Wrote simulated MP pipeline outputs to ", out_dir)
