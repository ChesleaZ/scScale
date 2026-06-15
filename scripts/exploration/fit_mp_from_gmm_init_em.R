library(ggplot2)

gmm_dir <- "outputs/exploration/gmm_spectrum"
out_dir <- "outputs/exploration/mp_from_gmm_init_em"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_iter <- 300
tol <- 1e-8
density_floor <- 1e-12

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
  if (sum(w) <= 1e-10) {
    return(mean(x))
  }
  sum(w * x) / sum(w)
}

weighted_var <- function(x, w, mu) {
  w[!is.finite(w)] <- 0
  if (sum(w) <= 1e-10 || !is.finite(mu)) {
    return(var(x))
  }
  sum(w * (x - mu)^2) / sum(w)
}

moment_start <- function(vals, weights) {
  mu <- weighted_mean(vals, weights)
  v <- weighted_var(vals, weights, mu)
  lambda <- v / max(mu^2, .Machine$double.eps)
  if (!is.finite(lambda)) {
    lambda <- 1
  }
  if (!is.finite(mu) || mu <= 0) {
    mu <- mean(vals[vals > 0])
  }
  c(lambda = min(max(lambda, 1e-4), 50), sigma2 = max(mu, .Machine$double.eps))
}

fit_mp_weighted <- function(vals, weights, start) {
  weights[!is.finite(weights)] <- 0
  if (sum(weights) <= 1e-10) {
    weights <- rep(1, length(vals))
  }
  start <- c(lambda = as.numeric(start["lambda"]), sigma2 = as.numeric(start["sigma2"]))
  if (any(!is.finite(start)) || start["lambda"] <= 0 || start["sigma2"] <= 0) {
    start <- moment_start(vals, weights)
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
    c(lambda = start["lambda"], sigma2 = start["sigma2"] * 0.5),
    c(lambda = start["lambda"], sigma2 = start["sigma2"] * 2),
    c(lambda = 0.25, sigma2 = weighted_mean(vals, weights)),
    c(lambda = 1, sigma2 = weighted_mean(vals, weights)),
    c(lambda = 3, sigma2 = weighted_mean(vals, weights))
  )
  starts[, "lambda"] <- pmin(pmax(starts[, "lambda"], 1e-4), 100)
  starts[, "sigma2"] <- pmax(starts[, "sigma2"], .Machine$double.eps)

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = log(starts[i, ]),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1200)
    )
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  out <- c(lambda = as.numeric(exp(best$par[1])), sigma2 = as.numeric(exp(best$par[2])), objective = best$value)
  if (any(!is.finite(out)) || out["lambda"] <= 0 || out["sigma2"] <= 0) {
    out <- c(lambda = as.numeric(start["lambda"]), sigma2 = as.numeric(start["sigma2"]), objective = NA_real_)
  }
  out
}

fit_two_mp_em <- function(vals, init_resp) {
  vals <- vals[is.finite(vals) & vals > 0]
  resp <- init_resp
  resp <- resp / rowSums(resp)
  starts <- list(
    moment_start(vals, resp[, 1]),
    moment_start(vals, resp[, 2])
  )

  loglik_trace <- numeric(max_iter)
  param_trace <- matrix(NA_real_, nrow = max_iter, ncol = 6)
  colnames(param_trace) <- c("pi1", "pi2", "lambda1", "lambda2", "sigma1", "sigma2")

  for (iter in seq_len(max_iter)) {
    fit1 <- fit_mp_weighted(vals, resp[, 1], starts[[1]])
    fit2 <- fit_mp_weighted(vals, resp[, 2], starts[[2]])
    fits <- list(fit1, fit2)

    ord <- order(c(fit1["sigma2"], fit2["sigma2"]))
    fits <- fits[ord]
    resp <- resp[, ord, drop = FALSE]
    starts <- lapply(fits, function(fit) c(lambda = as.numeric(fit["lambda"]), sigma2 = as.numeric(fit["sigma2"])))

    pi <- colMeans(resp)
    d1 <- pmax(mp_positive_density(vals, fits[[1]]["lambda"], fits[[1]]["sigma2"]), density_floor)
    d2 <- pmax(mp_positive_density(vals, fits[[2]]["lambda"], fits[[2]]["sigma2"]), density_floor)
    weighted_density <- cbind(pi[1] * d1, pi[2] * d2)
    total_density <- pmax(rowSums(weighted_density), density_floor)
    resp <- weighted_density / total_density
    resp[!is.finite(resp)] <- 0
    empty <- rowSums(resp) <= 0
    resp[empty, ] <- matrix(rep(c(0.5, 0.5), sum(empty)), ncol = 2, byrow = TRUE)
    resp <- resp / rowSums(resp)

    pi <- colMeans(resp)
    total_density <- pmax(pi[1] * d1 + pi[2] * d2, density_floor)
    loglik_trace[iter] <- sum(log(total_density))
    param_trace[iter, ] <- c(
      pi,
      fits[[1]]["lambda"],
      fits[[2]]["lambda"],
      fits[[1]]["sigma2"],
      fits[[2]]["sigma2"]
    )

    if (
      iter > 1 &&
        is.finite(loglik_trace[iter]) &&
        abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol * (1 + abs(loglik_trace[iter - 1]))
    ) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      param_trace <- param_trace[seq_len(iter), , drop = FALSE]
      break
    }
  }

  final <- param_trace[nrow(param_trace), ]
  lambda <- c(final["lambda1"], final["lambda2"])
  sigma2 <- c(final["sigma1"], final["sigma2"])
  weight <- c(final["pi1"], final["pi2"])

  list(
    resp = resp,
    parameters = data.frame(
      component = c("mp_low", "mp_high"),
      weight = weight,
      p_hat = weight * length(vals),
      lambda = lambda,
      sigma2 = sigma2,
      mp_lower = sigma2 * (1 - sqrt(lambda))^2,
      mp_upper = sigma2 * (1 + sqrt(lambda))^2,
      n_eigenvalues = length(vals),
      loglik = tail(loglik_trace, 1),
      n_iter = length(loglik_trace)
    ),
    loglik_trace = data.frame(iter = seq_along(loglik_trace), loglik = loglik_trace)
  )
}

make_diagnostics <- function(vals, fit) {
  params <- fit$parameters
  x_grid <- seq(max(min(vals), 1e-12), max(vals, params$mp_upper, na.rm = TRUE), length.out = 20000)
  curve_df <- rbind(
    data.frame(component = params$component[1], eigenvalue = x_grid, density = params$weight[1] * mp_positive_density(x_grid, params$lambda[1], params$sigma2[1])),
    data.frame(component = params$component[2], eigenvalue = x_grid, density = params$weight[2] * mp_positive_density(x_grid, params$lambda[2], params$sigma2[2]))
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

  qq_probs <- seq(0.0025, 0.9975, length.out = 600)
  qq_df <- data.frame(
    prob = qq_probs,
    empirical = as.numeric(quantile(vals, qq_probs, names = FALSE)),
    model = approx(mix_df$cdf, mix_df$eigenvalue, xout = qq_probs, ties = "ordered", rule = 2)$y
  )

  log_grid <- seq(log10(min(vals)), log10(max(vals)), length.out = 20000)
  raw_from_log <- 10^log_grid
  log_curve_df <- rbind(
    data.frame(component = params$component[1], log10_eigenvalue = log_grid, density = params$weight[1] * mp_positive_density(raw_from_log, params$lambda[1], params$sigma2[1]) * raw_from_log * log(10)),
    data.frame(component = params$component[2], log10_eigenvalue = log_grid, density = params$weight[2] * mp_positive_density(raw_from_log, params$lambda[2], params$sigma2[2]) * raw_from_log * log(10))
  )
  log_mix_df <- aggregate(density ~ log10_eigenvalue, log_curve_df, sum)

  write.csv(mix_df, file.path(out_dir, "mp_from_gmm_model_cdf.csv"), row.names = FALSE)
  write.csv(empirical_df, file.path(out_dir, "mp_from_gmm_empirical_vs_model_cdf.csv"), row.names = FALSE)
  write.csv(qq_df, file.path(out_dir, "mp_from_gmm_qq.csv"), row.names = FALSE)

  theme_set(theme_classic(base_size = 13))

  raw_density_plot <- ggplot(data.frame(eigenvalue = vals), aes(eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = curve_df, aes(eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = mix_df, aes(eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "MP EM initialized from GMM, outliers removed",
      subtitle = "Raw eigenvalue density; black curve is total MP mixture",
      x = "Covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  log_density_plot <- ggplot(data.frame(log10_eigenvalue = log10(vals)), aes(log10_eigenvalue)) +
    geom_histogram(aes(y = after_stat(density)), bins = 120, fill = "grey80", color = "white", linewidth = 0.1) +
    geom_line(data = log_curve_df, aes(log10_eigenvalue, density, color = component), linewidth = 0.8) +
    geom_line(data = log_mix_df, aes(log10_eigenvalue, density), color = "black", linewidth = 1.1) +
    labs(
      title = "MP EM initialized from GMM, log10 view",
      subtitle = "MP density transformed with the log10 Jacobian",
      x = "log10 covariance eigenvalue",
      y = "Density",
      color = "Component"
    )

  cdf_plot <- ggplot(empirical_df, aes(eigenvalue)) +
    geom_step(aes(y = empirical_cdf), color = "grey35", linewidth = 0.8) +
    geom_line(aes(y = model_cdf), color = "#1F78B4", linewidth = 0.9) +
    labs(
      title = "MP EM CDF, initialized from GMM",
      subtitle = "Outlier component removed before fitting",
      x = "Covariance eigenvalue",
      y = "CDF"
    )

  qq_lim <- range(c(qq_df$model, qq_df$empirical), finite = TRUE)
  qq_plot <- ggplot(qq_df, aes(model, empirical)) +
    geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(size = 1.35, alpha = 0.75, color = "#1F78B4") +
    coord_equal(xlim = qq_lim, ylim = qq_lim) +
    labs(
      title = "MP EM Q-Q, initialized from GMM",
      x = "Fitted MP mixture quantile",
      y = "Empirical quantile"
    )

  ggsave(file.path(out_dir, "mp_from_gmm_density_raw.png"), raw_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "mp_from_gmm_density_log10.png"), log_density_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "mp_from_gmm_cdf.png"), cdf_plot, width = 8.5, height = 5.25, dpi = 220)
  ggsave(file.path(out_dir, "mp_from_gmm_qq.png"), qq_plot, width = 6, height = 6, dpi = 220)
}

gmm_assign <- read.csv(file.path(gmm_dir, "log10_gmm_assignments.csv"))
gmm_params <- read.csv(file.path(gmm_dir, "log10_parameters.csv"))

outlier_col <- grep("^log_outlier$", names(gmm_assign), value = TRUE)
if (length(outlier_col) != 1) {
  stop("Expected log_outlier column in GMM assignments. Re-run fit_gmm_spectrum.R first.")
}

keep <- gmm_assign[[outlier_col]] < 0.5
bulk <- gmm_assign[keep, ]
vals <- bulk$eigenvalue
init_resp <- as.matrix(bulk[, c("log_g1", "log_g2")])
init_resp <- init_resp / rowSums(init_resp)

fit <- fit_two_mp_em(vals, init_resp)

assign_df <- data.frame(
  n_cells = bulk$n_cells,
  rank = bulk$rank,
  eigenvalue = vals,
  log10_eigenvalue = log10(vals),
  gmm_component = bulk$log_component,
  mp_low = fit$resp[, 1],
  mp_high = fit$resp[, 2]
)
assign_df$mp_component <- fit$parameters$component[max.col(fit$resp, ties.method = "first")]

initial_summary <- data.frame(
  component = c("gmm_g1", "gmm_g2", "removed_outlier"),
  p_hat = c(sum(init_resp[, 1]), sum(init_resp[, 2]), sum(!keep)),
  source = "log10_gmm_plus_outlier"
)

write.csv(gmm_params, file.path(out_dir, "gmm_initial_parameters.csv"), row.names = FALSE)
write.csv(initial_summary, file.path(out_dir, "gmm_initial_bulk_counts.csv"), row.names = FALSE)
write.csv(fit$parameters, file.path(out_dir, "mp_from_gmm_parameters.csv"), row.names = FALSE)
write.csv(assign_df, file.path(out_dir, "mp_from_gmm_assignments.csv"), row.names = FALSE)
write.csv(fit$loglik_trace, file.path(out_dir, "mp_from_gmm_loglik.csv"), row.names = FALSE)

make_diagnostics(vals, fit)

print(initial_summary)
print(fit$parameters)
message("Wrote MP EM initialized from GMM outputs to ", out_dir)
