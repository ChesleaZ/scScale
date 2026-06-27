# Fresh object-style interface for the scaling-law Gaussian spike fit.

scscale_fit <- function(
  x,
  input = c("counts", "normalized", "representation", "eigenvalues"),
  target = NULL,
  target_fit = NULL,
  n = NULL,
  p = NULL,
  target_depth = 1e4,
  center = TRUE,
  scale = TRUE,
  n_features = NULL,
  min_cells = 10,
  mp_max_iter = 50,
  mp_grid_n = 20000,
  mp_stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  mp_stop_tol = 1e-4,
  mp_min_iter = 2,
  r = NULL,
  fit_umi = input[1] == "counts",
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00),
  U_grid = NULL,
  umi_seed = 1,
  umi_replicates = 1,
  store_matrix = FALSE,
  use_irlba = TRUE
) {
  input <- match.arg(input)
  mp_stop_metric <- match.arg(mp_stop_metric)
  matched_call <- match.call()

  prepared <- scscale_prepare_matrix(
    x,
    input = input,
    target_depth = target_depth,
    center = center,
    scale = scale,
    n_features = n_features,
    min_cells = min_cells,
    n = n,
    p = p
  )

  spectrum <- scscale_spectrum(prepared$matrix, input = input, n = prepared$n, p = prepared$p)
  bulk <- scscale_fit_mp_iterative(
    spectrum$lambda,
    c_X = prepared$p / prepared$n,
    max_iter = mp_max_iter,
    grid_n = mp_grid_n,
    stop_metric = mp_stop_metric,
    stop_tol = mp_stop_tol,
    min_iter = mp_min_iter
  )
  spikes <- scscale_spike_table(spectrum$lambda, bulk, r = r)
  recoverability <- scscale_recoverability(spikes$q_X, c_X = prepared$p / prepared$n)
  recoverability_infinity <- scscale_theta_infinity(spikes$q_X)
  spike_keep <- if (is.null(r)) seq_len(nrow(spikes)) else seq_len(min(r, nrow(spikes)))

  object <- list(
    call = matched_call,
    input = input,
    n = prepared$n,
    p = prepared$p,
    c_X = prepared$p / prepared$n,
    c_ratio = prepared$p / prepared$n,
    normalization = prepared$normalization,
    features = prepared$features,
    spectrum = spectrum,
    bulk = bulk,
    spikes = spikes,
    theta_X = recoverability,
    theta = recoverability,
    theta_infinity = recoverability_infinity,
    r = length(spike_keep)
  )

  if (isTRUE(store_matrix)) {
    object$matrix <- prepared$matrix
  }

  if (!is.null(target)) {
    object$target <- scscale_target_mi(prepared$matrix, target, r = object$r, use_irlba = use_irlba)
  }

  if (!is.null(target_fit)) {
    object$mi <- scscale_mi(object, target_fit, r = object$r)
    object$theory_mi <- object$mi
    object$I_infinity <- object$mi$I_infinity
  }

  class(object) <- "scscale_fit"

  if (isTRUE(fit_umi)) {
    if (input != "counts") {
      stop("fit_umi = TRUE requires input = 'counts'.", call. = FALSE)
    }
    object$umi_scaling <- scscale_umi_scaling(
      x,
      U_grid = U_grid,
      sampling_rates = sampling_rates,
      r = object$r,
      target_depth = target_depth,
      n_features = n_features,
      min_cells = min_cells,
      seed = umi_seed,
      n_replicates = umi_replicates,
      theta_Y = if (is.null(target_fit)) NULL else target_fit$theta_X,
      reference_fit = object,
      empirical = FALSE,
      use_irlba = use_irlba,
      mp_max_iter = mp_max_iter,
      mp_grid_n = mp_grid_n
    )
  }

  object
}

scscale_prepare_matrix <- function(
  x,
  input,
  target_depth = 1e4,
  center = TRUE,
  scale = TRUE,
  n_features = NULL,
  min_cells = 10,
  n = NULL,
  p = NULL
) {
  if (input == "eigenvalues") {
    if (is.null(n) || is.null(p)) {
      stop("n and p are required when input = 'eigenvalues'.", call. = FALSE)
    }
    ev <- clean_eigenvalues(x)
    mat <- structure(matrix(numeric(0), nrow = p, ncol = n), eigenvalues = ev)
    return(list(
      matrix = mat,
      n = n,
      p = p,
      features = NULL,
      normalization = list(method = "none", input = "eigenvalues")
    ))
  }

  mat <- counts_matrix(x)
  if (is.null(rownames(mat))) rownames(mat) <- paste0("feature_", seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("cell_", seq_len(ncol(mat)))

  features <- rownames(mat)
  if (input == "counts" && !is.null(n_features)) {
    features <- select_hvgs(mat, n_features = n_features, min_cells = min_cells)
    mat <- mat[features, , drop = FALSE]
  }

  if (input == "counts") {
    mat <- scscale_normalize_counts(mat, target_depth = target_depth, center = center, scale = scale)
    norm <- list(
      method = "library_size_log1p_gene_scale",
      target_depth = target_depth,
      center = center,
      scale = scale
    )
  } else {
    mat <- as_dense_matrix(mat)
    if (isTRUE(center) || isTRUE(scale)) {
      mat <- scscale_standardize_rows(mat, center = center, scale = scale)
    }
    norm <- list(
      method = if (input == "representation") "representation_row_scale" else "row_scale",
      center = center,
      scale = scale
    )
  }

  list(
    matrix = mat,
    n = ncol(mat),
    p = nrow(mat),
    features = rownames(mat),
    normalization = norm
  )
}

scscale_normalize_counts <- function(counts, target_depth = 1e4, center = TRUE, scale = TRUE) {
  counts <- as_dense_matrix(counts)
  counts[counts < 0] <- 0
  library_size <- colSums(counts)
  library_size[!is.finite(library_size) | library_size <= 0] <- 1
  normalized <- sweep(counts, 2, library_size, "/") * target_depth
  normalized <- log1p(normalized)
  scscale_standardize_rows(normalized, center = center, scale = scale)
}

scscale_standardize_rows <- function(x, center = TRUE, scale = TRUE, scale_floor = 1e-12) {
  x <- as_dense_matrix(x)
  if (isTRUE(center)) {
    x <- x - rowMeans(x)
  }
  if (isTRUE(scale)) {
    row_sd <- apply(x, 1, stats::sd)
    row_sd[!is.finite(row_sd) | row_sd < scale_floor] <- 1
    x <- sweep(x, 1, row_sd, "/")
  }
  x[!is.finite(x)] <- 0
  x
}

scscale_spectrum <- function(x, input = "normalized", n = NULL, p = NULL) {
  ev_attr <- attr(x, "eigenvalues")
  if (!is.null(ev_attr)) {
    ev <- clean_eigenvalues(ev_attr)
    return(data.frame(
      rank = seq_along(ev),
      lambda = ev,
      eigenvalue = ev,
      log10_eigenvalue = log10(pmax(ev, .Machine$double.xmin)),
      n = n,
      p = p,
      c_X = p / n,
      source = input
    ))
  }

  x <- as_dense_matrix(x)
  n <- ncol(x)
  p <- nrow(x)
  if (n < 1L || p < 1L) {
    stop("x must have at least one row and one column.", call. = FALSE)
  }
  sv <- svd(t(x), nu = 0, nv = 0)$d
  ev <- sort((sv^2) / n, decreasing = TRUE)
  ev <- ev[is.finite(ev) & ev > 0]
  data.frame(
    rank = seq_along(ev),
    lambda = ev,
    eigenvalue = ev,
    log10_eigenvalue = log10(pmax(ev, .Machine$double.xmin)),
    n = n,
    p = p,
    c_X = p / n,
    source = input
  )
}

scscale_mp_edges <- function(c_X, tau2 = 1) {
  if (!is.finite(c_X) || !is.finite(tau2) || c_X <= 0 || tau2 <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(
    lower = tau2 * (1 - sqrt(c_X))^2,
    upper = tau2 * (1 + sqrt(c_X))^2
  )
}

scscale_mp_density <- function(x, c_X, tau2 = 1, positive_only = TRUE) {
  edge <- scscale_mp_edges(c_X, tau2)
  out <- numeric(length(x))
  inside <- x >= max(edge[["lower"]], 0) & x <= edge[["upper"]] & x > 0
  out[inside] <- sqrt((edge[["upper"]] - x[inside]) * (x[inside] - edge[["lower"]])) /
    (2 * pi * c_X * tau2 * x[inside])
  if (isTRUE(positive_only) && c_X > 1) out <- out * c_X
  out
}

scscale_mp_quantile <- function(c_X, prob = 0.5, tau2 = 1, grid_n = 20000) {
  edge <- scscale_mp_edges(c_X, tau2)
  if (!all(is.finite(edge))) return(NA_real_)
  x <- seq(max(edge[["lower"]], 0) + 1e-12, edge[["upper"]] - 1e-12, length.out = grid_n)
  density <- scscale_mp_density(x, c_X, tau2 = tau2, positive_only = TRUE)
  dx <- c(diff(x), utils::tail(diff(x), 1))
  cdf <- cumsum(density * dx)
  cdf <- cdf / max(cdf)
  stats::approx(cdf, x, xout = prob, ties = "ordered", rule = 2)$y
}

scscale_mp_cdf <- function(q, c_X, tau2 = 1, grid_n = 20000) {
  edge <- scscale_mp_edges(c_X, tau2)
  out <- rep(NA_real_, length(q))
  if (!all(is.finite(edge))) return(out)

  lower <- max(edge[["lower"]], 0) + 1e-12
  upper <- edge[["upper"]] - 1e-12
  x <- seq(lower, upper, length.out = grid_n)
  density <- scscale_mp_density(x, c_X, tau2 = tau2, positive_only = TRUE)
  dx <- c(diff(x), utils::tail(diff(x), 1))
  cdf <- cumsum(density * dx)
  cdf <- cdf / max(cdf)

  out[q <= lower] <- 0
  out[q >= upper] <- 1
  inside <- q > lower & q < upper
  out[inside] <- stats::approx(x, cdf, xout = q[inside], ties = "ordered", rule = 2)$y
  out
}

scscale_mp_goodness <- function(
  eigenvalues,
  c_X,
  tau2,
  grid_n = 20000,
  eps = 1e-12
) {
  ev <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0])
  m <- length(ev)
  if (m < 2L) {
    return(data.frame(
      qq_rmse_log = NA_real_,
      ks = NA_real_,
      neg_log_likelihood = NA_real_,
      support_fraction = NA_real_
    ))
  }

  edge <- scscale_mp_edges(c_X, tau2)
  prob <- (seq_len(m) - 0.5) / m
  mp_q <- vapply(
    prob,
    scscale_mp_quantile,
    numeric(1),
    c_X = c_X,
    tau2 = tau2,
    grid_n = grid_n
  )
  qq_rmse_log <- sqrt(mean((log(pmax(ev, eps)) - log(pmax(mp_q, eps)))^2))

  mp_cdf <- scscale_mp_cdf(ev, c_X = c_X, tau2 = tau2, grid_n = grid_n)
  ecdf_mid <- prob
  ks <- max(abs(ecdf_mid - mp_cdf), na.rm = TRUE)

  density <- scscale_mp_density(ev, c_X = c_X, tau2 = tau2, positive_only = TRUE)
  neg_log_likelihood <- mean(-log(pmax(density, eps)))
  support_fraction <- mean(ev >= max(edge[["lower"]], 0) & ev <= edge[["upper"]])

  data.frame(
    qq_rmse_log = qq_rmse_log,
    ks = ks,
    neg_log_likelihood = neg_log_likelihood,
    support_fraction = support_fraction
  )
}

scscale_fit_mp_iterative <- function(
  eigenvalues,
  c_X = NULL,
  c_ratio = NULL,
  max_iter = 50,
  grid_n = 20000,
  stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  stop_tol = 1e-4,
  min_iter = 2
) {
  stop_metric <- match.arg(stop_metric)
  if (is.null(c_X)) c_X <- c_ratio
  if (is.null(c_X)) stop("c_X is required.", call. = FALSE)
  ev <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
  if (length(ev) < 2L) {
    stop("Need at least two positive eigenvalues to fit the MP bulk.", call. = FALSE)
  }

  mp_median <- scscale_mp_quantile(c_X, prob = 0.5, tau2 = 1, grid_n = grid_n)
  bulk_index <- seq_along(ev)
  history <- vector("list", max_iter)
  best_metric <- Inf
  best_state <- NULL
  stop_reason <- "max_iter"

  for (iter in seq_len(max_iter)) {
    tau2 <- stats::median(ev[bulk_index]) / mp_median
    edge <- scscale_mp_edges(c_X, tau2)
    next_bulk <- which(ev <= edge[["upper"]] & ev > 0)
    spike_index <- which(ev > edge[["upper"]])
    gof <- scscale_mp_goodness(ev[next_bulk], c_X = c_X, tau2 = tau2, grid_n = grid_n)
    metric_value <- if (stop_metric == "none") NA_real_ else gof[[stop_metric]]
    history[[iter]] <- data.frame(
      iteration = iter,
      tau2 = tau2,
      c_X = c_X,
      mp_lower = edge[["lower"]],
      lambda_plus = edge[["upper"]],
      mp_upper = edge[["upper"]],
      n_bulk = length(next_bulk),
      n_spikes = length(spike_index),
      qq_rmse_log = gof$qq_rmse_log,
      ks = gof$ks,
      neg_log_likelihood = gof$neg_log_likelihood,
      support_fraction = gof$support_fraction,
      stop_metric = stop_metric,
      stop_metric_value = metric_value
    )

    meaningful_improvement <- FALSE
    if (stop_metric != "none" && is.finite(metric_value)) {
      meaningful_improvement <- is.null(best_state) || metric_value < best_metric - stop_tol
    }
    if (stop_metric != "none" && is.finite(metric_value) &&
        (is.null(best_state) || metric_value < best_metric)) {
      best_metric <- metric_value
      best_state <- list(
        tau2 = tau2,
        edge = edge,
        bulk_index = next_bulk,
        spike_index = spike_index,
        iteration = iter
      )
    }

    if (identical(next_bulk, bulk_index)) {
      stop_reason <- "bulk_stable"
      if (is.null(best_state)) {
        best_state <- list(
          tau2 = tau2,
          edge = edge,
          bulk_index = next_bulk,
          spike_index = spike_index,
          iteration = iter
        )
      }
      break
    }
    if (!length(next_bulk)) {
      stop_reason <- "empty_bulk"
      break
    }
    if (stop_metric != "none" && iter >= min_iter && is.finite(metric_value) && !meaningful_improvement &&
        !is.null(best_state) && metric_value > best_metric - stop_tol) {
      stop_reason <- "goodness_not_improving"
      break
    }
    bulk_index <- next_bulk
  }

  history <- do.call(rbind, history[seq_len(iter)])
  if (is.null(best_state)) {
    best_state <- list(
      tau2 = tau2,
      edge = edge,
      bulk_index = next_bulk,
      spike_index = spike_index,
      iteration = iter
    )
  }
  tau2 <- best_state$tau2
  edge <- best_state$edge
  bulk_index <- best_state$bulk_index
  spike_index <- best_state$spike_index
  list(
    method = "iterative_median_mp",
    tau2 = tau2,
    c_X = c_X,
    c_ratio = c_X,
    mp_median = mp_median,
    mp_lower = edge[["lower"]],
    lambda_plus = edge[["upper"]],
    mp_upper = edge[["upper"]],
    bulk_index = bulk_index,
    spike_index = spike_index,
    iterations = iter,
    selected_iteration = best_state$iteration,
    stop_metric = stop_metric,
    stop_reason = stop_reason,
    history = history
  )
}

scscale_spike_table <- function(eigenvalues, bulk, r = NULL) {
  ev <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
  lambda_tilde <- ev / bulk$tau2
  term <- lambda_tilde - (1 + bulk$c_X)
  disc <- term^2 - 4 * bulk$c_X
  d2 <- rep(0, length(ev))
  ok <- seq_along(ev) %in% bulk$spike_index & is.finite(disc) & disc >= 0
  d2[ok] <- pmax((term[ok] + sqrt(disc[ok])) / (2 * bulk$c_X), 0)
  q <- bulk$c_X * d2

  out <- data.frame(
    rank = seq_along(ev),
    lambda = ev,
    eigenvalue = ev,
    lambda_tilde = lambda_tilde,
    is_spike = seq_along(ev) %in% bulk$spike_index,
    d2_X = d2,
    d2 = d2,
    d = sqrt(d2),
    q_X = q,
    q = q,
    tau2 = bulk$tau2,
    c_X = bulk$c_X,
    c_ratio = bulk$c_X,
    lambda_plus = bulk$lambda_plus,
    mp_upper = bulk$mp_upper
  )

  if (!is.null(r)) out <- out[seq_len(min(r, nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

scscale_recoverability <- function(q, c_X = NULL, c_ratio = NULL) {
  if (is.null(c_X)) c_X <- c_ratio
  if (is.null(c_X)) {
    stop("c_X is required.", call. = FALSE)
  }
  theta <- rep(0, length(q))
  ok <- is.finite(q) & q > 0 & q^2 > c_X
  theta[ok] <- (q[ok]^2 - c_X) / (q[ok] * (q[ok] + 1))
  pmin(pmax(theta, 0), 1 - 1e-12)
}

scscale_theta_infinity <- function(q) {
  theta <- rep(0, length(q))
  ok <- is.finite(q) & q > 0
  theta[ok] <- q[ok] / (q[ok] + 1)
  pmin(pmax(theta, 0), 1 - 1e-12)
}

scscale_target_mi <- function(x, target, r = NULL, eps = 1e-12, use_irlba = TRUE) {
  r_use <- r %||% min(10L, ncol(x) - 1L, nrow(x) - 1L)
  empirical_mi(x, target, r = r_use, eps = eps, use_irlba = use_irlba)
}

scscale_theory_mi <- function(fit, target_fit, r = NULL, eps = 1e-12) {
  if (!inherits(fit, "scscale_fit") || !inherits(target_fit, "scscale_fit")) {
    stop("fit and target_fit must be scscale_fit objects.", call. = FALSE)
  }
  r_use <- r %||% min(length(fit$theta_X), length(target_fit$theta_X))
  theta_X <- fit$theta_X[seq_len(min(r_use, length(fit$theta_X)))]
  theta_Y <- target_fit$theta_X[seq_len(min(r_use, length(target_fit$theta_X)))]
  m <- min(length(theta_X), length(theta_Y))
  theta_double <- pmin(pmax(theta_X[seq_len(m)] * theta_Y[seq_len(m)], 0), 1 - eps)
  I_theory <- -0.5 * sum(log1p(-theta_double))
  list(
    I_theory = I_theory,
    mi = I_theory,
    theta_X = theta_X[seq_len(m)],
    theta_Y = theta_Y[seq_len(m)],
    theta_x = theta_X[seq_len(m)],
    theta_y = theta_Y[seq_len(m)],
    theta_double = theta_double,
    r = m
  )
}

scscale_i_infinity <- function(x, y = NULL, r = NULL, eps = 1e-12) {
  q_X <- if (inherits(x, "scscale_fit")) x$spikes$q_X else as.numeric(x)
  q_Y <- if (inherits(y, "scscale_fit")) y$spikes$q_X else as.numeric(y)
  if (is.null(q_Y)) {
    stop("Provide y as a scscale_fit object or numeric q_Y vector.", call. = FALSE)
  }

  r_use <- r %||% min(length(q_X), length(q_Y))
  m <- min(r_use, length(q_X), length(q_Y))
  q_X <- q_X[seq_len(m)]
  q_Y <- q_Y[seq_len(m)]
  theta_X_infinity <- scscale_theta_infinity(q_X)
  theta_Y_infinity <- scscale_theta_infinity(q_Y)
  theta_double_infinity <- pmin(
    pmax(theta_X_infinity * theta_Y_infinity, 0),
    1 - eps
  )

  list(
    I_infinity = -0.5 * sum(log1p(-theta_double_infinity)),
    mi = -0.5 * sum(log1p(-theta_double_infinity)),
    q_X = q_X,
    q_Y = q_Y,
    theta_X_infinity = theta_X_infinity,
    theta_Y_infinity = theta_Y_infinity,
    theta_double_infinity = theta_double_infinity,
    r = m
  )
}

scscale_mi <- function(
  fit,
  target_fit,
  n_grid = NULL,
  U_grid = NULL,
  sampling_rates = NULL,
  combine = TRUE,
  r = NULL,
  eps = 1e-12
) {
  if (!inherits(fit, "scscale_fit") || !inherits(target_fit, "scscale_fit")) {
    stop("fit and target_fit must be scscale_fit objects.", call. = FALSE)
  }

  r_use <- r %||% min(length(fit$theta_X), length(target_fit$theta_X))
  m <- min(r_use, nrow(fit$spikes), nrow(target_fit$spikes))
  theta_Y <- target_fit$theta_X[seq_len(m)]
  theta_Y_infinity <- target_fit$theta_infinity[seq_len(m)]

  base <- scscale_information_from_theta(fit$theta_X[seq_len(m)], theta_Y, eps = eps)
  infinite <- scscale_information_from_theta(
    fit$theta_infinity[seq_len(m)],
    theta_Y_infinity,
    eps = eps
  )

  out <- list(
    I = base,
    I_theory = base,
    mi = base,
    I_infinity = infinite,
    theta_X = fit$theta_X[seq_len(m)],
    theta_Y = theta_Y,
    theta_X_infinity = fit$theta_infinity[seq_len(m)],
    theta_Y_infinity = theta_Y_infinity,
    r = m
  )

  if (!is.null(n_grid)) {
    out$cell_scaling <- stats::aggregate(
      I_theory ~ n + c_X,
      data = scscale_cell_scaling(
        fit$spikes$d2_X[seq_len(m)],
        p = fit$p,
        n_grid = n_grid,
        theta_Y = theta_Y,
        eps = eps
      ),
      FUN = unique
    )
  }

  umi <- fit$umi_scaling
  if (!is.null(U_grid) || !is.null(sampling_rates)) {
    if (is.null(umi)) {
      stop("fit does not contain UMI scaling. Refit with fit_umi = TRUE or call scscale_umi_scaling().", call. = FALSE)
    }
    if (!is.null(sampling_rates)) {
      umi_rows <- umi$scaling[umi$scaling$sampling_rate %in% sampling_rates, , drop = FALSE]
    } else {
      umi_rows <- umi$scaling[umi$scaling$U %in% U_grid, , drop = FALSE]
    }
    umi_summary <- unique(umi_rows[, c("sampling_rate", "U", "I_theory", "I_infinity")])
    out$umi_scaling <- umi_summary

    if (isTRUE(combine) && !is.null(n_grid) && nrow(umi_rows)) {
      rate_values <- sort(unique(umi_rows$sampling_rate))
      combined <- do.call(rbind, lapply(rate_values, function(rate) {
        rate_rows <- umi_rows[umi_rows$sampling_rate == rate, , drop = FALSE]
        q_df <- unique(rate_rows[, c("rank", "q_X")])
        q_df <- q_df[order(q_df$rank), , drop = FALSE]
        d2_X <- q_df$q_X[seq_len(min(m, nrow(q_df)))] / fit$c_X
        curve <- stats::aggregate(
          I_theory ~ n + c_X,
          data = scscale_cell_scaling(d2_X, p = fit$p, n_grid = n_grid, theta_Y = theta_Y, eps = eps),
          FUN = unique
        )
        curve$sampling_rate <- rate
        curve
      }))
      out$grid <- combined
    }
  }

  class(out) <- "scscale_mi"
  out
}

print.scscale_fit <- function(x, ...) {
  cat("scscale_fit\n")
  cat("  input: ", x$input, "\n", sep = "")
  cat("  dimensions: p=", x$p, ", n=", x$n, ", c_X=", signif(x$c_X, 4), "\n", sep = "")
  cat("  MP: tau2=", signif(x$bulk$tau2, 4), ", lambda_plus=", signif(x$bulk$lambda_plus, 4),
      ", iterations=", x$bulk$iterations, "\n", sep = "")
  cat("  spikes: ", sum(x$spikes$is_spike), "\n", sep = "")
  if (!is.null(x$theory_mi)) cat("  theory MI: ", signif(x$theory_mi$mi, 4), "\n", sep = "")
  if (!is.null(x$I_infinity)) cat("  I_infinity: ", signif(x$I_infinity, 4), "\n", sep = "")
  if (!is.null(x$target)) cat("  target MI: ", signif(x$target$mi, 4), "\n", sep = "")
  invisible(x)
}
