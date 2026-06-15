library(ggplot2)

in_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
out_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
clip_at <- 0.5
n_bins <- 180

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

fit_mp <- function(eigvals) {
  vals <- eigvals[eigvals > 1e-10 & eigvals <= clip_at]
  probs <- seq(0.02, 0.98, by = 0.02)
  empirical_q <- as.numeric(quantile(vals, probs, names = FALSE))

  objective <- function(par) {
    gamma <- exp(par[1])
    sigma2 <- exp(par[2])
    fitted_q <- mp_quantiles(gamma, sigma2, probs)
    mean((log(empirical_q) - log(fitted_q))^2)
  }

  starts <- expand.grid(
    gamma = c(0.2, 0.5, 1, 2, 5),
    sigma2 = c(median(vals), mean(vals), quantile(vals, 0.75))
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
    gamma = gamma,
    sigma2 = sigma2,
    lower = sigma2 * (1 - sqrt(gamma))^2,
    upper = sigma2 * (1 + sqrt(gamma))^2,
    objective = best$value,
    n_fit = length(vals)
  )
}

eig_df <- read.csv(in_file)
all_n <- max(eig_df$n_cells)
all_df <- subset(eig_df, n_cells == all_n)
fit <- fit_mp(all_df$eigenvalue)

x <- seq(max(fit$lower, 1e-12), min(fit$upper, clip_at), length.out = 1000)
curve_df <- data.frame(
  eigenvalue = x,
  density = mp_density(x, fit$gamma, fit$sigma2)
)

all_df$eigenvalue_clipped <- pmin(all_df$eigenvalue, clip_at)
n_gt_clip <- sum(all_df$eigenvalue > clip_at)

write.csv(
  fit,
  file.path(out_dir, "jurkat_all_cells_clip_0_5_mp_fit.csv"),
  row.names = FALSE
)

p <- ggplot(all_df, aes(eigenvalue_clipped)) +
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
    linewidth = 1
  ) +
  geom_vline(xintercept = clip_at, color = "#C43C39", linetype = "dashed") +
  annotate(
    "text",
    x = clip_at * 0.96,
    y = Inf,
    label = paste0(">", clip_at, ": ", n_gt_clip),
    hjust = 1,
    vjust = 1.4,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "All-cell Jurkat eigenvalue distribution",
    subtitle = paste0("Eigenvalues clipped at ", clip_at, "; free MP fit over positive clipped range"),
    x = "Covariance eigenvalue, clipped",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "jurkat_all_cells_clip_0_5_mp_fit.png"),
  p,
  width = 8,
  height = 5,
  dpi = 220
)

print(fit)
cat("n_cells =", all_n, "\n")
cat("n_eigenvalues_gt_0.5 =", n_gt_clip, "\n")
