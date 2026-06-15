library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/two_mp_free_lambda_em"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cells_use <- NULL
max_iter <- 300
tol <- 1e-7
density_floor <- 1e-12
spike_tail_quantile <- 0.95

mp_density <- function(x, lambda, sigma2) {
  lower <- sigma2 * (1 - sqrt(lambda))^2
  upper <- sigma2 * (1 + sqrt(lambda))^2
  out <- numeric(length(x))
  inside <- x >= lower & x <= upper & x > 0
  out[inside] <- sqrt((upper - x[inside]) * (x[inside] - lower)) /
    (2 * pi * lambda * sigma2 * x[inside])
  out
}

mp_positive_mass <- function(lambda) {
  ifelse(lambda > 1, 1 / lambda, 1)
}

mp_positive_density <- function(x, lambda, sigma2) {
  mp_density(x, lambda, sigma2) / mp_positive_mass(lambda)
}

fit_mp_weighted <- function(vals, weights) {
  positive_vals <- vals[vals > 0]
  q <- quantile(positive_vals, c(0.05, 0.95), names = FALSE)

  objective <- function(par) {
    lambda <- exp(par[1])
    sigma2 <- exp(par[2])
    d <- pmax(mp_positive_density(vals, lambda, sigma2), density_floor)
    -sum(weights * log(d))
  }

  starts <- expand.grid(
    lambda = c(0.02, 0.05, 0.1, 0.25, 0.5, 0.9, 1.5, 3, 8),
    sigma2 = unique(as.numeric(c(mean(positive_vals), median(positive_vals), q)))
  )

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(c(starts$lambda[i], starts$sigma2[i])),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1000)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  c(lambda = exp(best$par[1]), sigma2 = exp(best$par[2]), objective = best$value)
}

fit_em <- function(vals) {
  vals <- vals[is.finite(vals) & vals > 0]
  p <- length(vals)
  spike_lower <- as.numeric(quantile(vals, spike_tail_quantile))
  spike_density_value <- 1 / (max(vals) - spike_lower)

  resp <- cbind(
    as.numeric(vals <= quantile(vals, 0.55)),
    as.numeric(vals > quantile(vals, 0.30) & vals < spike_lower),
    as.numeric(vals >= spike_lower)
  )
  resp <- resp / rowSums(resp)

  loglik_trace <- numeric(max_iter)
  params <- matrix(NA_real_, nrow = max_iter, ncol = 7)
  colnames(params) <- c("pi1", "pi2", "pi3", "lambda1", "lambda2", "sigma1", "sigma2")

  for (iter in seq_len(max_iter)) {
    pi <- colMeans(resp)

    fit1 <- fit_mp_weighted(vals, resp[, 1])
    fit2 <- fit_mp_weighted(vals, resp[, 2])
    ord <- order(c(fit1["sigma2"], fit2["sigma2"]))
    fits <- list(fit1, fit2)[ord]
    resp[, 1:2] <- resp[, ord]
    pi <- colMeans(resp)

    d1 <- pmax(mp_positive_density(vals, fits[[1]]["lambda"], fits[[1]]["sigma2"]), density_floor)
    d2 <- pmax(mp_positive_density(vals, fits[[2]]["lambda"], fits[[2]]["sigma2"]), density_floor)
    d3 <- ifelse(vals >= spike_lower, spike_density_value, density_floor)

    weighted_density <- cbind(pi[1] * d1, pi[2] * d2, pi[3] * d3)
    total_density <- rowSums(weighted_density)
    resp <- weighted_density / total_density

    loglik_trace[iter] <- sum(log(total_density))
    params[iter, ] <- c(
      colMeans(resp),
      fits[[1]]["lambda"],
      fits[[2]]["lambda"],
      fits[[1]]["sigma2"],
      fits[[2]]["sigma2"]
    )

    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      params <- params[seq_len(iter), , drop = FALSE]
      break
    }
  }

  final <- params[nrow(params), ]
  sigma2 <- c(final["sigma1"], final["sigma2"])
  lambda <- c(final["lambda1"], final["lambda2"])

  list(
    resp = resp,
    parameters = data.frame(
      component = c("noise_mp", "intrinsic_mp", "spike_outlier"),
      weight = c(final["pi1"], final["pi2"], final["pi3"]),
      p_hat = p * c(final["pi1"], final["pi2"], final["pi3"]),
      lambda = c(lambda, NA),
      sigma2 = c(sigma2, NA),
      mp_lower = c(sigma2 * (1 - sqrt(lambda))^2, NA),
      mp_upper = c(sigma2 * (1 + sqrt(lambda))^2, NA),
      spike_lower = spike_lower,
      n_eigenvalues = p,
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace)
  )
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
n_cells <- if (is.null(n_cells_use)) max(eig_df$n_cells) else n_cells_use
eigvals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells], decreasing = TRUE)
fit <- fit_em(eigvals)

assign_df <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  noise_mp = fit$resp[, 1],
  intrinsic_mp = fit$resp[, 2],
  spike_outlier = fit$resp[, 3]
)
assign_df$component <- c("noise_mp", "intrinsic_mp", "spike_outlier")[max.col(fit$resp, ties.method = "first")]

write.csv(fit$parameters, file.path(out_dir, "two_mp_free_lambda_em_parameters.csv"), row.names = FALSE)
write.csv(assign_df, file.path(out_dir, "two_mp_free_lambda_em_assignments.csv"), row.names = FALSE)
write.csv(fit$loglik_trace, file.path(out_dir, "two_mp_free_lambda_em_loglik.csv"), row.names = FALSE)

x_max <- as.numeric(quantile(eigvals, 0.995))
params <- fit$parameters
spike_lower <- params$spike_lower[1]
spike_density_value <- 1 / (max(eigvals) - spike_lower)
model_x_max <- max(eigvals, params$mp_upper[1:2], na.rm = TRUE)
x_model <- seq(max(min(eigvals), 1e-12), model_x_max, length.out = 20000)

curve_df <- rbind(
  data.frame(component = "noise_mp", eigenvalue = x_model, density = params$weight[1] * mp_positive_density(x_model, params$lambda[1], params$sigma2[1])),
  data.frame(component = "intrinsic_mp", eigenvalue = x_model, density = params$weight[2] * mp_positive_density(x_model, params$lambda[2], params$sigma2[2])),
  data.frame(component = "spike_outlier", eigenvalue = x_model, density = ifelse(x_model >= spike_lower & x_model <= max(eigvals), params$weight[3] * spike_density_value, 0))
)
mix_df <- aggregate(density ~ eigenvalue, curve_df, sum)
dx <- c(diff(mix_df$eigenvalue), tail(diff(mix_df$eigenvalue), 1))
mix_df$cdf <- cumsum(mix_df$density * dx)
mix_df$cdf <- mix_df$cdf / max(mix_df$cdf)

empirical_df <- data.frame(
  eigenvalue = sort(eigvals),
  empirical_cdf = seq_along(eigvals) / length(eigvals)
)
empirical_df$model_cdf <- approx(mix_df$eigenvalue, mix_df$cdf, xout = empirical_df$eigenvalue, rule = 2)$y

qq_probs <- seq(0.005, 0.995, length.out = 500)
qq_df <- data.frame(
  prob = qq_probs,
  empirical = as.numeric(quantile(eigvals, qq_probs, names = FALSE)),
  model = approx(mix_df$cdf, mix_df$eigenvalue, xout = qq_probs, ties = "ordered", rule = 2)$y
)

write.csv(mix_df, file.path(out_dir, "two_mp_free_lambda_em_model_cdf.csv"), row.names = FALSE)
write.csv(empirical_df, file.path(out_dir, "two_mp_free_lambda_em_empirical_vs_model_cdf.csv"), row.names = FALSE)
write.csv(qq_df, file.path(out_dir, "two_mp_free_lambda_em_qq.csv"), row.names = FALSE)

density_plot_df <- data.frame(eigenvalue = eigvals[eigvals <= x_max])
density_norm <- mean(eigvals <= x_max)
density_curve_df <- subset(curve_df, eigenvalue <= x_max)
density_curve_df$density <- density_curve_df$density / density_norm
density_mix_df <- subset(mix_df, eigenvalue <= x_max)
density_mix_df$density <- density_mix_df$density / density_norm

theme_set(theme_classic(base_size = 13))

density_plot <- ggplot(density_plot_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 140, fill = "grey78", color = "white", linewidth = 0.1) +
  geom_line(data = density_curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
  geom_line(data = density_mix_df, aes(eigenvalue, density), linewidth = 1.1, color = "black") +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(
    title = "Two-MP free-lambda EM fit",
    subtitle = paste0("Density is renormalized over displayed eigenvalues <= ", signif(x_max, 3)),
    x = "Covariance eigenvalue",
    y = "Density",
    color = "Component"
  )

scree_plot <- ggplot(assign_df, aes(rank, eigenvalue, color = component)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Two-MP free-lambda EM assignments",
    x = "Eigenvalue rank",
    y = "Covariance eigenvalue",
    color = "MAP component"
  )

loglik_plot <- ggplot(fit$loglik_trace, aes(iter, loglik)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(title = "EM log likelihood", x = "Iteration", y = "Log likelihood")

cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
  geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
  geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(
    title = "Two-MP free-lambda EM CDF check",
    subtitle = "Grey: empirical CDF; blue: fitted mixture CDF",
    x = "Covariance eigenvalue",
    y = "CDF"
  )

qq_plot <- ggplot(qq_df, aes(model, empirical)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
  geom_point(size = 1.4, alpha = 0.75, color = "#1F78B4") +
  coord_equal(xlim = c(0, x_max), ylim = c(0, x_max)) +
  labs(
    title = "Two-MP free-lambda EM Q-Q check",
    x = "Fitted mixture quantile",
    y = "Empirical eigenvalue quantile"
  )

ggsave(file.path(out_dir, "two_mp_free_lambda_em_density_fit.png"), density_plot, width = 9, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_free_lambda_em_scree_assignments.png"), scree_plot, width = 8, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_free_lambda_em_loglik.png"), loglik_plot, width = 7, height = 4.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_free_lambda_em_cdf.png"), cdf_plot, width = 8, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "two_mp_free_lambda_em_qq.png"), qq_plot, width = 6, height = 6, dpi = 220)

print(fit$parameters)
message("Wrote two-MP free-lambda EM outputs to ", out_dir)
