library(ggplot2)

in_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
out_dir <- "outputs/exploration/gaussian_init_mp_em"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cells_use <- NULL
spike_tail_quantile <- 0.98
max_iter_gaussian <- 300
max_iter_mp <- 300
tol <- 1e-7
density_floor <- 1e-12

normal_density <- function(x, mean, sd) {
  if (!is.finite(mean) || !is.finite(sd) || sd <= 0) {
    return(rep(0, length(x)))
  }
  pmax(dnorm(x, mean = mean, sd = sd), density_floor)
}

mp_density <- function(x, lambda, sigma2) {
  if (!is.finite(lambda) || !is.finite(sigma2) || lambda <= 0 || sigma2 <= 0) {
    return(rep(0, length(x)))
  }
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

weighted_mean <- function(x, w) {
  w[!is.finite(w)] <- 0
  if (sum(w) <= 1e-8) {
    return(mean(x))
  }
  sum(w * x) / sum(w)
}

weighted_var <- function(x, w, mu) {
  w[!is.finite(w)] <- 0
  if (sum(w) <= 1e-8 || !is.finite(mu)) {
    return(var(x))
  }
  sum(w * (x - mu)^2) / sum(w)
}

fit_gaussian_init <- function(vals) {
  vals <- vals[is.finite(vals) & vals > 0]
  p <- length(vals)
  spike_lower <- as.numeric(quantile(vals, spike_tail_quantile))
  spike_density_value <- 1 / (max(vals) - spike_lower)

  q <- quantile(vals, c(0.50, 0.90), names = FALSE)
  resp <- cbind(
    as.numeric(vals <= q[1]),
    as.numeric(vals > q[1] & vals < spike_lower),
    as.numeric(vals >= spike_lower)
  )
  resp <- resp / rowSums(resp)

  loglik_trace <- numeric(max_iter_gaussian)
  params_trace <- matrix(NA_real_, nrow = max_iter_gaussian, ncol = 9)
  colnames(params_trace) <- c("pi1", "pi2", "pi3", "mean1", "mean2", "sd1", "sd2", "spike_lower", "spike_density")

  for (iter in seq_len(max_iter_gaussian)) {
    pi <- colMeans(resp)
    mu <- c(weighted_mean(vals, resp[, 1]), weighted_mean(vals, resp[, 2]))
    sd <- sqrt(pmax(c(
      weighted_var(vals, resp[, 1], mu[1]),
      weighted_var(vals, resp[, 2], mu[2])
    ), 1e-10))

    ord <- order(mu)
    mu <- mu[ord]
    sd <- sd[ord]
    resp[, 1:2] <- resp[, ord]
    pi <- colMeans(resp)

    d1 <- normal_density(vals, mu[1], sd[1])
    d2 <- normal_density(vals, mu[2], sd[2])
    d3 <- ifelse(vals >= spike_lower, spike_density_value, density_floor)

    weighted_density <- cbind(pi[1] * d1, pi[2] * d2, pi[3] * d3)
    weighted_density[!is.finite(weighted_density)] <- 0
    total_density <- pmax(rowSums(weighted_density), density_floor)
    resp <- weighted_density / total_density
    resp[!is.finite(resp)] <- 0
    empty <- rowSums(resp) <= 0
    resp[empty, ] <- matrix(rep(c(1 / 3, 1 / 3, 1 / 3), sum(empty)), ncol = 3, byrow = TRUE)
    resp <- resp / rowSums(resp)

    loglik_trace[iter] <- sum(log(total_density))
    params_trace[iter, ] <- c(colMeans(resp), mu, sd, spike_lower, spike_density_value)

    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      params_trace <- params_trace[seq_len(iter), , drop = FALSE]
      break
    }
  }

  final <- as.numeric(params_trace[nrow(params_trace), ])
  names(final) <- colnames(params_trace)
  list(
    resp = resp,
    parameters = data.frame(
      component = c("gaussian_low", "gaussian_high", "spike_outlier"),
      weight = final[c("pi1", "pi2", "pi3")],
      p_hat = p * final[c("pi1", "pi2", "pi3")],
      mean = c(final["mean1"], final["mean2"], NA),
      variance = c(final["sd1"]^2, final["sd2"]^2, NA),
      sd = c(final["sd1"], final["sd2"], NA),
      spike_lower = final["spike_lower"],
      n_eigenvalues = p,
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace),
    spike_lower = final["spike_lower"],
    spike_density_value = final["spike_density"]
  )
}

moment_mp_start <- function(vals, weights) {
  weights[!is.finite(weights)] <- 0
  if (sum(weights) <= 1e-8) {
    weights <- rep(1, length(vals))
  }
  mu <- weighted_mean(vals, weights)
  v <- weighted_var(vals, weights, mu)
  lambda <- v / max(mu^2, .Machine$double.eps)
  if (!is.finite(lambda)) {
    lambda <- 1
  }
  if (!is.finite(mu) || mu <= 0) {
    mu <- mean(vals[vals > 0])
  }
  c(lambda = min(max(lambda, 1e-4), 20), sigma2 = max(mu, .Machine$double.eps))
}

fit_mp_weighted <- function(vals, weights, start) {
  weights[!is.finite(weights)] <- 0
  if (sum(weights) <= 1e-8) {
    weights <- rep(1, length(vals))
  }
  start <- c(lambda = as.numeric(start[1]), sigma2 = as.numeric(start[2]))
  if (any(!is.finite(start)) || start["lambda"] <= 0 || start["sigma2"] <= 0) {
    start <- moment_mp_start(vals, weights)
  }
  objective <- function(par) {
    lambda <- exp(par[1])
    sigma2 <- exp(par[2])
    d <- pmax(mp_positive_density(vals, lambda, sigma2), density_floor)
    -sum(weights * log(d))
  }

  starts <- rbind(
    start,
    c(lambda = start["lambda"] * 0.5, sigma2 = start["sigma2"]),
    c(lambda = start["lambda"] * 2, sigma2 = start["sigma2"]),
    c(lambda = start["lambda"], sigma2 = start["sigma2"] * 0.75),
    c(lambda = start["lambda"], sigma2 = start["sigma2"] * 1.25)
  )
  starts[, "lambda"] <- pmin(pmax(starts[, "lambda"], 1e-4), 50)
  starts[, "sigma2"] <- pmax(starts[, "sigma2"], .Machine$double.eps)

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(starts[i, ]),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1000)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  out <- c(lambda = as.numeric(exp(best$par[1])), sigma2 = as.numeric(exp(best$par[2])), objective = best$value)
  if (any(!is.finite(out)) || out["lambda"] <= 0 || out["sigma2"] <= 0) {
    out <- c(lambda = as.numeric(start["lambda"]), sigma2 = as.numeric(start["sigma2"]), objective = NA_real_)
  }
  out
}

fit_mp_refinement <- function(vals, gaussian_fit) {
  vals <- vals[is.finite(vals) & vals > 0]
  p <- length(vals)
  spike_lower <- gaussian_fit$spike_lower
  spike_density_value <- gaussian_fit$spike_density_value
  resp <- gaussian_fit$resp
  spike_resp <- resp[, 3]
  non_spike_mass <- pmax(1 - spike_resp, 0)

  starts <- list(
    moment_mp_start(vals, resp[, 1]),
    moment_mp_start(vals, resp[, 2])
  )

  loglik_trace <- numeric(max_iter_mp)
  params_trace <- matrix(NA_real_, nrow = max_iter_mp, ncol = 9)
  colnames(params_trace) <- c("pi1", "pi2", "pi3", "lambda1", "lambda2", "sigma1", "sigma2", "spike_lower", "spike_density")

  for (iter in seq_len(max_iter_mp)) {
    pi <- colMeans(resp)

    fit1 <- fit_mp_weighted(vals, resp[, 1], starts[[1]])
    fit2 <- fit_mp_weighted(vals, resp[, 2], starts[[2]])
    fits <- list(fit1, fit2)

    ord <- order(c(fit1["sigma2"], fit2["sigma2"]))
    fits <- fits[ord]
    resp[, 1:2] <- resp[, ord]
    starts <- lapply(fits, function(fit) c(lambda = as.numeric(fit["lambda"]), sigma2 = as.numeric(fit["sigma2"])))
    pi <- colMeans(resp)

    d1 <- pmax(mp_positive_density(vals, fits[[1]]["lambda"], fits[[1]]["sigma2"]), density_floor)
    d2 <- pmax(mp_positive_density(vals, fits[[2]]["lambda"], fits[[2]]["sigma2"]), density_floor)
    d3 <- ifelse(vals >= spike_lower, spike_density_value, density_floor)

    mp_pi <- colMeans(resp[, 1:2, drop = FALSE])
    mp_weighted_density <- cbind(mp_pi[1] * d1, mp_pi[2] * d2)
    mp_weighted_density[!is.finite(mp_weighted_density)] <- 0
    mp_total_density <- pmax(rowSums(mp_weighted_density), density_floor)
    mp_resp <- mp_weighted_density / mp_total_density
    mp_resp[!is.finite(mp_resp)] <- 0
    empty <- rowSums(mp_resp) <= 0
    mp_resp[empty, ] <- matrix(rep(c(1 / 2, 1 / 2), sum(empty)), ncol = 2, byrow = TRUE)
    mp_resp <- mp_resp / rowSums(mp_resp)

    resp[, 1:2] <- non_spike_mass * mp_resp
    resp[, 3] <- spike_resp

    pi <- colMeans(resp)
    pi[!is.finite(pi)] <- 0
    pi <- pi / sum(pi)
    total_density <- pmax(pi[1] * d1 + pi[2] * d2 + pi[3] * d3, density_floor)
    loglik_trace[iter] <- sum(log(total_density))
    params_trace[iter, ] <- c(
      colMeans(resp),
      fits[[1]]["lambda"],
      fits[[2]]["lambda"],
      fits[[1]]["sigma2"],
      fits[[2]]["sigma2"],
      spike_lower,
      spike_density_value
    )

    if (
      iter > 1 &&
        is.finite(loglik_trace[iter]) &&
        is.finite(loglik_trace[iter - 1]) &&
        abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))
    ) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      params_trace <- params_trace[seq_len(iter), , drop = FALSE]
      break
    }
  }

  final <- as.numeric(params_trace[nrow(params_trace), ])
  names(final) <- colnames(params_trace)
  lambda <- final[c("lambda1", "lambda2")]
  sigma2 <- final[c("sigma1", "sigma2")]

  list(
    resp = resp,
    parameters = data.frame(
      component = c("noise_mp", "intrinsic_mp", "spike_outlier"),
      weight = final[c("pi1", "pi2", "pi3")],
      p_hat = p * final[c("pi1", "pi2", "pi3")],
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

model_diagnostics <- function(vals, parameters, prefix) {
  x_max <- as.numeric(quantile(vals, 0.995))
  spike_lower <- parameters$spike_lower[1]
  spike_density_value <- 1 / (max(vals) - spike_lower)
  model_x_max <- max(vals, parameters$mp_upper[1:2], na.rm = TRUE)
  x_model <- seq(max(min(vals), 1e-12), model_x_max, length.out = 20000)

  curve_df <- rbind(
    data.frame(component = "noise_mp", eigenvalue = x_model, density = parameters$weight[1] * mp_positive_density(x_model, parameters$lambda[1], parameters$sigma2[1])),
    data.frame(component = "intrinsic_mp", eigenvalue = x_model, density = parameters$weight[2] * mp_positive_density(x_model, parameters$lambda[2], parameters$sigma2[2])),
    data.frame(component = "spike_outlier", eigenvalue = x_model, density = ifelse(x_model >= spike_lower & x_model <= max(vals), parameters$weight[3] * spike_density_value, 0))
  )
  mix_df <- aggregate(density ~ eigenvalue, curve_df, sum)
  dx <- c(diff(mix_df$eigenvalue), tail(diff(mix_df$eigenvalue), 1))
  mix_df$cdf <- cumsum(mix_df$density * dx)
  mix_df$cdf <- mix_df$cdf / max(mix_df$cdf)

  empirical_df <- data.frame(
    eigenvalue = sort(vals),
    empirical_cdf = seq_along(vals) / length(vals)
  )
  empirical_df$model_cdf <- approx(mix_df$eigenvalue, mix_df$cdf, xout = empirical_df$eigenvalue, rule = 2)$y

  qq_probs <- seq(0.005, 0.995, length.out = 500)
  qq_df <- data.frame(
    prob = qq_probs,
    empirical = as.numeric(quantile(vals, qq_probs, names = FALSE)),
    model = approx(mix_df$cdf, mix_df$eigenvalue, xout = qq_probs, ties = "ordered", rule = 2)$y
  )

  write.csv(mix_df, file.path(out_dir, paste0(prefix, "_model_cdf.csv")), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, paste0(prefix, "_empirical_vs_model_cdf.csv")), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, paste0(prefix, "_qq.csv")), row.names = FALSE)

  density_plot_df <- data.frame(eigenvalue = vals[vals <= x_max])
  density_norm <- mean(vals <= x_max)
  density_curve_df <- subset(curve_df, eigenvalue <= x_max)
  density_curve_df$density <- density_curve_df$density / density_norm
  density_mix_df <- subset(mix_df, eigenvalue <= x_max)
  density_mix_df$density <- density_mix_df$density / density_norm

  density_plot <- ggplot(density_plot_df, aes(eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 140, fill = "grey78", color = "white", linewidth = 0.1) +
    geom_line(data = density_curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = density_mix_df, aes(eigenvalue, density), linewidth = 1.1, color = "black") +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(
      title = "Gaussian-initialized MP EM density",
      subtitle = paste0("Displayed density renormalized for eigenvalues <= ", signif(x_max, 3)),
      x = "Covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(
      title = "Gaussian-initialized MP EM CDF",
      subtitle = "Grey: empirical CDF; blue: fitted mixture CDF",
      x = "Covariance eigenvalue",
      y = "CDF"
    )

  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.4, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = c(0, x_max), ylim = c(0, x_max)) +
    labs(
      title = "Gaussian-initialized MP EM Q-Q",
      x = "Fitted mixture quantile",
      y = "Empirical eigenvalue quantile"
    )

  ggsave(file.path(out_dir, paste0(prefix, "_density_fit.png")), density_plot, width = 9, height = 5.5, dpi = 220)
  ggsave(file.path(out_dir, paste0(prefix, "_cdf.png")), cdf_plot, width = 8, height = 5.5, dpi = 220)
  ggsave(file.path(out_dir, paste0(prefix, "_qq.png")), qq_plot, width = 6, height = 6, dpi = 220)
}

eig_df <- read.csv(file.path(in_dir, "jurkat_eigenvalues_by_cells.csv"))
n_cells <- if (is.null(n_cells_use)) max(eig_df$n_cells) else n_cells_use
eigvals <- sort(eig_df$eigenvalue[eig_df$n_cells == n_cells], decreasing = TRUE)

gaussian_fit <- fit_gaussian_init(eigvals)
mp_fit <- fit_mp_refinement(eigvals, gaussian_fit)

gaussian_assign_df <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  gaussian_low = gaussian_fit$resp[, 1],
  gaussian_high = gaussian_fit$resp[, 2],
  spike_outlier = gaussian_fit$resp[, 3]
)
gaussian_assign_df$component <- c("gaussian_low", "gaussian_high", "spike_outlier")[max.col(gaussian_fit$resp, ties.method = "first")]

mp_assign_df <- data.frame(
  n_cells = n_cells,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  noise_mp = mp_fit$resp[, 1],
  intrinsic_mp = mp_fit$resp[, 2],
  spike_outlier = mp_fit$resp[, 3]
)
mp_assign_df$component <- c("noise_mp", "intrinsic_mp", "spike_outlier")[max.col(mp_fit$resp, ties.method = "first")]

write.csv(gaussian_fit$parameters, file.path(out_dir, "gaussian_init_parameters.csv"), row.names = FALSE)
write.csv(gaussian_assign_df, file.path(out_dir, "gaussian_init_assignments.csv"), row.names = FALSE)
write.csv(gaussian_fit$loglik_trace, file.path(out_dir, "gaussian_init_loglik.csv"), row.names = FALSE)
write.csv(mp_fit$parameters, file.path(out_dir, "mp_refined_parameters.csv"), row.names = FALSE)
write.csv(mp_assign_df, file.path(out_dir, "mp_refined_assignments.csv"), row.names = FALSE)
write.csv(mp_fit$loglik_trace, file.path(out_dir, "mp_refined_loglik.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 13))

model_diagnostics(eigvals, mp_fit$parameters, "mp_refined")

assignment_plot <- ggplot(mp_assign_df, aes(rank, eigenvalue, color = component)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Gaussian-initialized MP EM assignments",
    x = "Eigenvalue rank",
    y = "Covariance eigenvalue",
    color = "MAP component"
  )

loglik_plot <- ggplot(mp_fit$loglik_trace, aes(iter, loglik)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(title = "MP refinement EM log likelihood", x = "Iteration", y = "Log likelihood")

ggsave(file.path(out_dir, "mp_refined_scree_assignments.png"), assignment_plot, width = 8, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "mp_refined_loglik.png"), loglik_plot, width = 7, height = 4.5, dpi = 220)

print(gaussian_fit$parameters)
print(mp_fit$parameters)
message("Wrote Gaussian initialization plus MP EM outputs to ", out_dir)
