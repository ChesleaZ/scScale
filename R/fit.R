# Fresh object-style interface for the scaling-law Gaussian spike fit.

scscale_fit <- function(
  x,
  input = c("counts", "normalized", "representation", "eigenvalues"),
  target = NULL,
  target_fit = NULL,
  n = NULL,
  p = NULL,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  n_features = NULL,
  min_cells = 10,
  mp_max_iter = 300,
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
  count_transform <- match.arg(count_transform)
  mp_stop_metric <- match.arg(mp_stop_metric)
  matched_call <- match.call()

  prepared <- scscale_prepare_matrix(
    x,
    input = input,
    target_depth = target_depth,
    count_transform = count_transform,
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
      count_transform = count_transform,
      center = center,
      scale = scale,
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
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  n_features = NULL,
  min_cells = 10,
  n = NULL,
  p = NULL
) {
  count_transform <- match.arg(count_transform)

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
    mat <- scscale_normalize_counts(
      mat,
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale
    )
    norm <- list(
      method = paste0(count_transform, "_row_center"),
      target_depth = target_depth,
      count_transform = count_transform,
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

scscale_normalize_counts <- function(
  counts,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE
) {
  count_transform <- match.arg(count_transform)
  counts <- as_dense_matrix(counts)
  counts[counts < 0] <- 0
  normalized <- switch(
    count_transform,
    log1p = log1p(counts),
    log1p_cpm = {
      library_size <- colSums(counts)
      library_size[!is.finite(library_size) | library_size <= 0] <- 1
      log1p(sweep(counts, 2, library_size, "/") * target_depth)
    },
    pearson_residual = pearson_residuals(counts)
  )
  scscale_standardize_rows(normalized, center = center, scale = scale)
}

scscale_standardize_rows <- function(x, center = TRUE, scale = FALSE, scale_floor = 1e-12) {
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

.scscale_mp_median_cache <- new.env(parent = emptyenv())

scscale_mp_cdf_integrated <- function(q, c_X, tau2 = 1) {
  edge <- scscale_mp_edges(c_X, tau2)
  if (!all(is.finite(edge))) return(NA_real_)
  lower <- max(edge[["lower"]], 0)
  upper <- edge[["upper"]]
  if (q <= lower) return(0)
  if (q >= upper) return(1)

  stats::integrate(
    scscale_mp_density,
    lower = lower,
    upper = q,
    c_X = c_X,
    tau2 = tau2,
    positive_only = TRUE,
    subdivisions = 200L,
    rel.tol = 1e-8
  )$value
}

scscale_mp_quantile_grid <- function(c_X, prob = 0.5, tau2 = 1, grid_n = 20000) {
  edge <- scscale_mp_edges(c_X, tau2)
  if (!all(is.finite(edge))) return(rep(NA_real_, length(prob)))
  lower <- max(edge[["lower"]], 0)
  upper <- edge[["upper"]]
  dx <- (upper - lower) / grid_n
  x <- lower + (seq_len(grid_n) - 0.5) * dx
  density <- scscale_mp_density(x, c_X, tau2 = tau2, positive_only = TRUE)
  cdf <- cumsum(density * dx)
  cdf <- cdf / max(cdf)
  stats::approx(cdf, x, xout = prob, ties = "ordered", rule = 2)$y
}

scscale_mp_quantile <- function(c_X, prob = 0.5, tau2 = 1, grid_n = 20000) {
  if (length(prob) != 1L) {
    return(scscale_mp_quantile_grid(c_X = c_X, prob = prob, tau2 = tau2, grid_n = grid_n))
  }

  edge <- scscale_mp_edges(c_X, tau2)
  if (!all(is.finite(edge))) return(NA_real_)
  lower <- max(edge[["lower"]], 0)
  upper <- edge[["upper"]]
  if (prob <= 0) return(lower)
  if (prob >= 1) return(upper)

  stats::uniroot(
    function(x) scscale_mp_cdf_integrated(x, c_X = c_X, tau2 = tau2) - prob,
    interval = c(lower, upper),
    tol = .Machine$double.eps^0.25
  )$root
}

scscale_mp_median_sample <- function(
  c_X,
  tau2 = 1,
  p_sim = 1200,
  n_reps = 3,
  seed = 1
) {
  if (!is.finite(c_X) || c_X <= 0 || !is.finite(tau2) || tau2 <= 0) {
    return(NA_real_)
  }
  p_sim <- max(20L, as.integer(p_sim))
  n_sim <- max(20L, as.integer(round(p_sim / c_X)))
  p_sim <- max(1L, as.integer(round(c_X * n_sim)))
  n_reps <- max(1L, as.integer(n_reps))
  cache_key <- paste(signif(c_X, 12), signif(tau2, 12), p_sim, n_sim, n_reps, seed, sep = ":")
  if (exists(cache_key, envir = .scscale_mp_median_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .scscale_mp_median_cache, inherits = FALSE))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  medians <- numeric(n_reps)
  for (rep_i in seq_len(n_reps)) {
    z <- matrix(stats::rnorm(p_sim * n_sim), nrow = p_sim, ncol = n_sim)
    sv <- svd(t(z), nu = 0, nv = 0)$d
    ev <- (sv^2) / n_sim
    ev <- ev[is.finite(ev) & ev > 0]
    medians[rep_i] <- stats::median(ev)
  }
  out <- tau2 * stats::median(medians)
  assign(cache_key, out, envir = .scscale_mp_median_cache)
  out
}

scscale_mp_cdf <- function(q, c_X, tau2 = 1, grid_n = 20000) {
  edge <- scscale_mp_edges(c_X, tau2)
  out <- rep(NA_real_, length(q))
  if (!all(is.finite(edge))) return(out)

  lower <- max(edge[["lower"]], 0)
  upper <- edge[["upper"]]
  dx <- (upper - lower) / grid_n
  x <- lower + (seq_len(grid_n) - 0.5) * dx
  density <- scscale_mp_density(x, c_X, tau2 = tau2, positive_only = TRUE)
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
  mp_q <- scscale_mp_quantile_grid(c_X = c_X, prob = prob, tau2 = tau2, grid_n = grid_n)
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
  max_iter = 300,
  grid_n = 20000,
  stop_metric = c("qq_rmse_log", "ks", "neg_log_likelihood", "none"),
  stop_tol = 1e-4,
  min_iter = 2,
  median_method = c("sample", "integrate", "grid"),
  mp_sample_p = 1200,
  mp_sample_reps = 3,
  mp_sample_seed = 1,
  gap_z_threshold = 5,
  quiet_z_threshold = 5,
  quiet_run = 20,
  background_skip = 25,
  background_window = 500,
  central_quantiles = c(0.10, 0.70)
) {
  stop_metric <- match.arg(stop_metric)
  median_method <- match.arg(median_method)
  if (is.null(c_X)) c_X <- c_ratio
  if (is.null(c_X)) stop("c_X is required.", call. = FALSE)
  ev <- sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
  if (length(ev) < 2L) {
    stop("Need at least two positive eigenvalues to fit the MP bulk.", call. = FALSE)
  }

  mp_median <- switch(
    median_method,
    sample = scscale_mp_median_sample(
      c_X,
      tau2 = 1,
      p_sim = mp_sample_p,
      n_reps = mp_sample_reps,
      seed = mp_sample_seed
    ),
    integrate = scscale_mp_quantile(c_X, prob = 0.5, tau2 = 1, grid_n = grid_n),
    grid = scscale_mp_quantile_grid(c_X, prob = 0.5, tau2 = 1, grid_n = grid_n)
  )

  k_max <- min(as.integer(max_iter), length(ev) - 20L)
  if (k_max < 1L) {
    tau2 <- stats::median(ev) / mp_median
    edge <- scscale_mp_edges(c_X, tau2)
    gof <- scscale_mp_goodness(ev, c_X = c_X, tau2 = tau2, grid_n = grid_n)
    history <- data.frame(
      iteration = 1L,
      K = 0L,
      tau2 = tau2,
      c_X = c_X,
      mp_lower = edge[["lower"]],
      lambda_plus = edge[["upper"]],
      mp_upper = edge[["upper"]],
      n_bulk = length(ev),
      n_spikes = 0L,
      log_gap = NA_real_,
      background_gap_median = NA_real_,
      background_gap_mad = NA_real_,
      gap_z = NA_real_,
      qq_rmse_log = gof$qq_rmse_log,
      ks = gof$ks,
      neg_log_likelihood = gof$neg_log_likelihood,
      support_fraction = gof$support_fraction,
      stop_metric = stop_metric,
      stop_metric_value = if (stop_metric == "none") NA_real_ else gof[[stop_metric]]
    )
    return(list(
      method = "median_mp_mad_gap",
      tau2 = tau2,
      c_X = c_X,
      c_ratio = c_X,
      mp_median = mp_median,
      mp_lower = edge[["lower"]],
      lambda_plus = edge[["upper"]],
      mp_upper = edge[["upper"]],
      bulk_index = seq_along(ev),
      spike_index = integer(0),
      iterations = 1L,
      selected_iteration = 1L,
      stop_metric = stop_metric,
      stop_reason = "too_few_candidates",
      history = history
    ))
  }

  log_ev <- log(pmax(ev, 1e-12))
  log_gaps <- log_ev[-length(log_ev)] - log_ev[-1L]
  history <- vector("list", k_max)
  for (K in seq_len(k_max)) {
    bulk_index <- seq.int(K + 1L, length(ev))
    tau2 <- stats::median(ev[bulk_index]) / mp_median
    edge <- scscale_mp_edges(c_X, tau2)
    gof <- scscale_mp_goodness(ev[bulk_index], c_X = c_X, tau2 = tau2, grid_n = grid_n)

    prob <- (seq_along(bulk_index) - 0.5) / length(bulk_index)
    mp_q <- scscale_mp_quantile_grid(c_X = c_X, prob = prob, tau2 = tau2, grid_n = grid_n)
    bulk_inc <- sort(ev[bulk_index], decreasing = FALSE)
    keep <- prob >= central_quantiles[1] & prob <= central_quantiles[2]
    central_loss <- stats::median((log(pmax(bulk_inc[keep], 1e-12)) - log(pmax(mp_q[keep], 1e-12)))^2)

    bg_start <- min(length(log_gaps), K + background_skip)
    bg_end <- min(length(log_gaps), K + background_skip + background_window - 1L)
    background_gaps <- log_gaps[bg_start:bg_end]
    background_median <- stats::median(background_gaps, na.rm = TRUE)
    background_mad <- stats::median(abs(background_gaps - background_median), na.rm = TRUE)
    gap_z <- (log_gaps[K] - background_median) / (1.4826 * background_mad + 1e-12)

    metric_value <- if (stop_metric == "none") NA_real_ else gof[[stop_metric]]
    history[[K]] <- data.frame(
      iteration = K,
      K = K,
      tau2 = tau2,
      c_X = c_X,
      mp_lower = edge[["lower"]],
      lambda_plus = edge[["upper"]],
      mp_upper = edge[["upper"]],
      n_bulk = length(bulk_index),
      n_spikes = K,
      log_gap = log_gaps[K],
      background_gap_median = background_median,
      background_gap_mad = background_mad,
      gap_z = gap_z,
      central_qq_loss = central_loss,
      qq_rmse_log = gof$qq_rmse_log,
      ks = gof$ks,
      neg_log_likelihood = gof$neg_log_likelihood,
      support_fraction = gof$support_fraction,
      stop_metric = stop_metric,
      stop_metric_value = metric_value
    )
  }

  history <- do.call(rbind, history)
  quiet_after <- vapply(seq_len(nrow(history)), function(i) {
    j <- seq.int(i + 1L, min(nrow(history), i + quiet_run))
    length(j) > 0 && all(history$gap_z[j] < quiet_z_threshold, na.rm = TRUE)
  }, logical(1))
  strong <- is.finite(history$gap_z) & history$gap_z >= gap_z_threshold
  candidate <- history[strong & quiet_after, , drop = FALSE]
  if (nrow(candidate) > 0L) {
    selected <- candidate[which.max(candidate$K), , drop = FALSE]
    stop_reason <- "last_strong_gap_before_quiet_region"
  } else if (any(strong)) {
    candidate <- history[strong, , drop = FALSE]
    selected <- candidate[which.max(candidate$K), , drop = FALSE]
    stop_reason <- "last_strong_gap"
  } else {
    selected <- history[which.max(history$gap_z), , drop = FALSE]
    stop_reason <- "largest_gap_z_fallback"
  }

  tau2 <- selected$tau2
  edge <- c(lower = selected$mp_lower, upper = selected$lambda_plus)
  spike_index <- seq_len(selected$K)
  bulk_index <- seq.int(selected$K + 1L, length(ev))
  list(
    method = "median_mp_mad_gap",
    tau2 = tau2,
    c_X = c_X,
    c_ratio = c_X,
    mp_median = mp_median,
    mp_lower = edge[["lower"]],
    lambda_plus = edge[["upper"]],
    mp_upper = edge[["upper"]],
    bulk_index = bulk_index,
    spike_index = spike_index,
    iterations = k_max,
    selected_iteration = selected$iteration,
    stop_metric = stop_metric,
    stop_reason = stop_reason,
    history = history,
    gap_z_threshold = gap_z_threshold,
    quiet_z_threshold = quiet_z_threshold,
    quiet_run = quiet_run,
    background_skip = background_skip,
    background_window = background_window,
    central_quantiles = central_quantiles
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

scscale_mi <- function(
  fit,
  target_fit,
  n_grid = NULL,
  U_grid = NULL,
  sampling_rates = NULL,
  combine = TRUE,
  empirical = TRUE,
  store_empirical_subspaces = FALSE,
  use_irlba = TRUE,
  r = NULL,
  P = NULL,
  eps = 1e-12
) {
  if (!inherits(fit, "scscale_fit") || !inherits(target_fit, "scscale_fit")) {
    stop("fit and target_fit must be scscale_fit objects.", call. = FALSE)
  }

  r_use <- r %||% min(length(fit$theta_X), length(target_fit$theta_X))
  m <- min(r_use, nrow(fit$spikes), nrow(target_fit$spikes))
  theta_Y <- target_fit$theta_X[seq_len(m)]
  theta_Y_infinity <- target_fit$theta_infinity[seq_len(m)]

  base_details <- scscale_low_rank_mi(fit$theta_X[seq_len(m)], theta_Y, P = P, eps = eps)
  base <- base_details$mi
  theta_double <- base_details$gamma
  infinite_details <- scscale_low_rank_mi(
    fit$theta_infinity[seq_len(m)],
    theta_Y_infinity,
    P = P,
    eps = eps
  )
  infinite <- infinite_details$mi
  theta_double_infinity <- infinite_details$gamma

  out <- list(
    I = base,
    I_theory = base,
    mi = base,
    I_infinity = infinite,
    theta_X = fit$theta_X[seq_len(m)],
    theta_Y = theta_Y,
    theta_X_infinity = fit$theta_infinity[seq_len(m)],
    theta_Y_infinity = theta_Y_infinity,
    r = m,
    intermediate = list(
      q_X = fit$spikes$q_X[seq_len(m)],
      q_Y = target_fit$spikes$q_X[seq_len(m)],
      d2_X = fit$spikes$d2_X[seq_len(m)],
      d2_Y = target_fit$spikes$d2_X[seq_len(m)],
      theta_double = theta_double,
      theta_double_infinity = theta_double_infinity,
      sigma = base_details$sigma,
      sigma_infinity = infinite_details$sigma,
      P = if (is.null(P)) NULL else base_details$P,
      c_X = fit$c_X,
      c_Y = target_fit$c_X,
      lambda_plus_X = fit$bulk$lambda_plus,
      lambda_plus_Y = target_fit$bulk$lambda_plus,
      tau2_X = fit$bulk$tau2,
      tau2_Y = target_fit$bulk$tau2
    )
  )

  if (isTRUE(empirical)) {
    has_matrices <- !is.null(fit$matrix) && !is.null(target_fit$matrix)
    if (has_matrices) {
      z_X <- right_singular_vectors(fit$matrix, r = m, use_irlba = use_irlba)
      z_Y <- right_singular_vectors(target_fit$matrix, r = m, use_irlba = use_irlba)
      empirical_overlap <- subspace_overlap_mi(z_X, z_Y, eps = eps)
      out$I_empirical <- empirical_overlap$mi
      out$empirical <- list(
        available = TRUE,
        I_empirical = empirical_overlap$mi,
        mi = empirical_overlap$mi,
        gamma = empirical_overlap$gamma,
        r_eff = empirical_overlap$r_eff
      )
      if (isTRUE(store_empirical_subspaces)) {
        out$empirical$z_X <- z_X
        out$empirical$z_Y <- z_Y
      }
    } else {
      out$I_empirical <- NA_real_
      out$empirical <- list(
        available = FALSE,
        I_empirical = NA_real_,
        mi = NA_real_,
        gamma = numeric(0),
        r_eff = 0L,
        reason = "Empirical MI requires both fit objects to contain stored normalized matrices. Refit with store_matrix = TRUE."
      )
    }
  }

  if (!is.null(n_grid)) {
    out$cell_scaling <- stats::aggregate(
      I_theory ~ n + c_X,
      data = scscale_cell_scaling(
        fit$spikes$q_X[seq_len(m)],
        p = fit$p,
        n_grid = n_grid,
        theta_Y = theta_Y,
        P = P,
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
        q_X <- q_df$q_X[seq_len(min(m, nrow(q_df)))]
        curve <- stats::aggregate(
          I_theory ~ n + c_X,
          data = scscale_cell_scaling(q_X, p = fit$p, n_grid = n_grid, theta_Y = theta_Y, P = P, eps = eps),
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
