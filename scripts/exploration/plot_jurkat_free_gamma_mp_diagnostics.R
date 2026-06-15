library(ggplot2)

eig_file <- "outputs/exploration/jurkat_eigenvalues_by_cells/jurkat_eigenvalues_by_cells.csv"
fit_file <- "outputs/exploration/jurkat_mp_by_ncells/jurkat_mp_by_ncells_parameters.csv"
out_dir <- "outputs/exploration/jurkat_mp_by_ncells"

grid_n <- 50000
qq_probs <- seq(0.01, 0.99, by = 0.005)

mp_density <- function(x, gamma, sigma2) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * gamma * sigma2 * x[inside])
  out
}

mp_grid <- function(gamma, sigma2, grid_n = 50000) {
  lower <- sigma2 * (1 - sqrt(gamma))^2
  upper <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(max(lower, 0) + 1e-12, upper - 1e-12, length.out = grid_n)
  d <- mp_density(x, gamma, sigma2)
  dx <- c(diff(x), tail(diff(x), 1))
  cdf <- cumsum(d * dx)
  cdf <- cdf / max(cdf)
  data.frame(eigenvalue = x, density = d, cdf = cdf)
}

eig_df <- read.csv(eig_file)
fit_df <- read.csv(fit_file)
fit_df <- fit_df[fit_df$fit_type == "free_gamma", ]
n_values <- sort(unique(fit_df$n_cells))

cdf_rows <- list()
qq_rows <- list()
curve_rows <- list()

for (n_cells in n_values) {
  vals <- eig_df$eigenvalue[eig_df$n_cells == n_cells]
  vals <- sort(vals[is.finite(vals) & vals > 1e-10])
  fit <- fit_df[fit_df$n_cells == n_cells, ]
  grid <- mp_grid(fit$gamma, fit$sigma2, grid_n = grid_n)
  grid$n_cells <- n_cells
  curve_rows[[length(curve_rows) + 1]] <- grid

  cdf_rows[[length(cdf_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    eigenvalue = vals,
    empirical_cdf = seq_along(vals) / length(vals),
    model_cdf = approx(grid$eigenvalue, grid$cdf, xout = vals, rule = 2)$y
  )

  qq_rows[[length(qq_rows) + 1]] <- data.frame(
    n_cells = n_cells,
    prob = qq_probs,
    empirical = as.numeric(quantile(vals, qq_probs, names = FALSE)),
    model = approx(grid$cdf, grid$eigenvalue, xout = qq_probs, ties = "ordered", rule = 2)$y
  )
}

cdf_df <- do.call(rbind, cdf_rows)
qq_df <- do.call(rbind, qq_rows)
curve_df <- do.call(rbind, curve_rows)

label_levels <- paste0("n = ", format(n_values, big.mark = ","))
cdf_df$n_cells_label <- factor(cdf_df$n_cells, levels = n_values, labels = label_levels)
qq_df$n_cells_label <- factor(qq_df$n_cells, levels = n_values, labels = label_levels)
curve_df$n_cells_label <- factor(curve_df$n_cells, levels = n_values, labels = label_levels)

write.csv(cdf_df, file.path(out_dir, "jurkat_free_gamma_mp_cdf.csv"), row.names = FALSE)
write.csv(qq_df, file.path(out_dir, "jurkat_free_gamma_mp_qq.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

cdf_plot <- ggplot(cdf_df, aes(eigenvalue)) +
  geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.65) +
  geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.85) +
  scale_x_log10() +
  facet_wrap(~ n_cells_label, ncol = 2) +
  labs(
    title = "Jurkat free-gamma MP fit CDF across cell counts",
    subtitle = "Positive eigenvalues only; grey: empirical CDF, blue: fitted MP CDF",
    x = "Covariance eigenvalue, log10 scale",
    y = "CDF"
  )

ggsave(file.path(out_dir, "jurkat_free_gamma_mp_cdf.png"), cdf_plot, width = 10, height = 8, dpi = 220)

qq_plot <- ggplot(qq_df, aes(model, empirical)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
  geom_point(size = 1.1, alpha = 0.75, color = "#1F78B4") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ n_cells_label, scales = "free", ncol = 2) +
  labs(
    title = "Jurkat free-gamma MP fit Q-Q across cell counts",
    subtitle = "Positive eigenvalue quantiles; both axes log10",
    x = "Fitted MP quantile",
    y = "Empirical quantile"
  )

ggsave(file.path(out_dir, "jurkat_free_gamma_mp_qq.png"), qq_plot, width = 10, height = 8, dpi = 220)

fit_label <- fit_df[, c("n_cells", "gamma_theory", "gamma", "sigma2", "mp_lower", "mp_upper", "objective")]
print(fit_label)
message("Wrote free-gamma MP diagnostics to ", out_dir)
