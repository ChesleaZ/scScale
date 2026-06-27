# Scaling curves from fitted Gaussian spike strengths.

scscale_downsample_counts <- function(counts, U, seed = NULL) {
  counts <- as_dense_matrix(counts)
  counts[counts < 0] <- 0
  if (!is.null(seed)) set.seed(seed)

  library_size <- colSums(counts)
  prob <- rep(1, ncol(counts))
  ok <- is.finite(library_size) & library_size > U & library_size > 0
  prob[ok] <- U / library_size[ok]

  out <- counts
  for (j in seq_len(ncol(out))) {
    if (prob[j] < 1) {
      out[, j] <- stats::rbinom(nrow(out), size = out[, j], prob = prob[j])
    }
  }
  rownames(out) <- rownames(counts)
  colnames(out) <- colnames(counts)
  out
}

scscale_downsample_counts_fraction <- function(counts, fraction, seed = NULL) {
  counts <- as_dense_matrix(counts)
  counts[counts < 0] <- 0
  if (!is.finite(fraction) || fraction <= 0 || fraction > 1) {
    stop("fraction must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  out <- matrix(
    stats::rbinom(length(counts), size = as.vector(counts), prob = fraction),
    nrow = nrow(counts),
    ncol = ncol(counts),
    dimnames = dimnames(counts)
  )
  out
}

scscale_information_from_theta <- function(theta_X, theta_Y, eps = 1e-12) {
  m <- min(length(theta_X), length(theta_Y))
  if (m < 1L) return(0)
  theta_double <- pmin(pmax(theta_X[seq_len(m)] * theta_Y[seq_len(m)], 0), 1 - eps)
  -0.5 * sum(log1p(-theta_double))
}

scscale_umi_scaling <- function(
  counts,
  U_grid = NULL,
  sampling_rates = NULL,
  r = 10,
  target_depth = 1e4,
  n_features = NULL,
  min_cells = 10,
  seed = 1,
  n_replicates = 1,
  theta_Y = NULL,
  reference_fit = NULL,
  empirical = TRUE,
  use_irlba = TRUE,
  mp_max_iter = 50,
  mp_grid_n = 20000
) {
  if (is.null(U_grid) && is.null(sampling_rates)) {
    stop("Provide U_grid or sampling_rates.", call. = FALSE)
  }
  if (!is.null(sampling_rates)) {
    if (length(sampling_rates) < 2L) {
      stop("sampling_rates must contain at least two values.", call. = FALSE)
    }
    if (any(!is.finite(sampling_rates) | sampling_rates <= 0 | sampling_rates > 1)) {
      stop("sampling_rates must be in (0, 1].", call. = FALSE)
    }
  } else if (length(U_grid) < 2L) {
    stop("U_grid must contain at least two UMI depths.", call. = FALSE)
  }
  n_replicates <- as.integer(n_replicates)
  if (!is.finite(n_replicates) || n_replicates < 1L) {
    stop("n_replicates must be a positive integer.", call. = FALSE)
  }
  counts <- counts_matrix(counts)
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature_", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

  features <- rownames(counts)
  if (!is.null(n_features)) {
    features <- select_hvgs(counts, n_features = n_features, min_cells = min_cells)
  }
  counts <- counts[features, , drop = FALSE]
  total_umi <- sum(counts)
  if (!is.null(sampling_rates)) {
    U_grid <- sampling_rates * total_umi
  }

  if (is.null(reference_fit)) {
    reference_fit <- scscale_fit(
      counts,
      input = "counts",
      target_depth = target_depth,
      r = r,
      mp_max_iter = mp_max_iter,
      mp_grid_n = mp_grid_n,
      fit_umi = FALSE,
      store_matrix = TRUE
    )
  } else if (!inherits(reference_fit, "scscale_fit")) {
    stop("reference_fit must be a scscale_fit object.", call. = FALSE)
  }

  X_ref <- reference_fit$matrix
  if (is.null(X_ref)) {
    prepared_ref <- scscale_prepare_matrix(
      counts,
      input = "counts",
      target_depth = target_depth,
      center = TRUE,
      scale = TRUE
    )
    X_ref <- prepared_ref$matrix
  }

  r_eff <- min(as.integer(r), nrow(X_ref), ncol(X_ref))
  if (r_eff < 1L) stop("No reference spike directions are available.", call. = FALSE)
  u_ref <- svd(X_ref, nu = r_eff, nv = 0)$u[, seq_len(r_eff), drop = FALSE]
  z_ref <- NULL
  if (isTRUE(empirical)) {
    z_ref <- right_singular_vectors(X_ref, r = r_eff, drop_first = FALSE, use_irlba = use_irlba)
  }

  if (is.null(theta_Y)) {
    theta_Y <- reference_fit$theta_infinity[seq_len(min(r_eff, length(reference_fit$theta_infinity)))]
  }

  rows <- vector("list", length(U_grid) * n_replicates)
  row_i <- 0L
  for (i in seq_along(U_grid)) {
    U <- U_grid[i]
    sampling_rate <- if (is.null(sampling_rates)) NA_real_ else sampling_rates[i]
    for (replicate in seq_len(n_replicates)) {
      row_i <- row_i + 1L
      replicate_seed <- seed + (i - 1L) * n_replicates + replicate
      if (is.null(sampling_rates)) {
        counts_U <- scscale_downsample_counts(counts, U = U, seed = replicate_seed)
      } else {
        counts_U <- scscale_downsample_counts_fraction(counts, fraction = sampling_rate, seed = replicate_seed)
      }
      X_U <- scscale_normalize_counts(counts_U, target_depth = target_depth, center = TRUE, scale = TRUE)
      fit_U <- scscale_fit(
        X_U,
        input = "normalized",
        center = FALSE,
        scale = FALSE,
        r = r_eff,
        mp_max_iter = mp_max_iter,
        mp_grid_n = mp_grid_n,
        fit_umi = FALSE
      )

      projected <- crossprod(u_ref, X_U)
      v <- rowSums(projected^2) / ncol(X_U)
      delta <- pmax(v - fit_U$bulk$tau2, 0)
      d2_X <- delta / fit_U$bulk$tau2
      q_X <- fit_U$c_X * d2_X
      theta_X <- scscale_recoverability(q_X, c_X = fit_U$c_X)
      theta_X_infinity <- scscale_theta_infinity(q_X)
      I_theory <- scscale_information_from_theta(theta_X, theta_Y)
      I_infinity <- scscale_information_from_theta(theta_X_infinity, theta_Y)
      empirical_overlap <- NULL
      gamma_empirical <- rep(NA_real_, r_eff)
      I_empirical <- NA_real_
      if (isTRUE(empirical) && !is.null(z_ref)) {
        z_U <- right_singular_vectors(X_U, r = r_eff, drop_first = FALSE, use_irlba = use_irlba)
        if (!is.null(z_U)) {
          empirical_overlap <- subspace_overlap_mi(z_U, z_ref)
          I_empirical <- empirical_overlap$mi
          gamma_empirical[seq_len(min(r_eff, length(empirical_overlap$gamma)))] <-
            empirical_overlap$gamma[seq_len(min(r_eff, length(empirical_overlap$gamma)))]
        }
      }

      rows[[row_i]] <- data.frame(
        U = U,
        replicate = replicate,
        rank = seq_len(r_eff),
        v = v,
        tau2 = fit_U$bulk$tau2,
        d2_X = d2_X,
        q_X = q_X,
        theta_X = theta_X,
        theta_X_infinity = theta_X_infinity,
        I_theory = I_theory,
        I_infinity = I_infinity,
        I_empirical = I_empirical,
        gamma_empirical = gamma_empirical,
        c_X = fit_U$c_X,
        n = fit_U$n,
        p = fit_U$p,
        sampling_rate = sampling_rate,
        total_umi_expected = U,
        total_umi_observed = sum(counts_U)
      )
    }
  }

  scaling <- do.call(rbind, rows)
  q_fit <- do.call(rbind, lapply(split(scaling, scaling$rank), function(df) {
    model <- stats::lm(q_X ~ U, data = df)
    data.frame(
      rank = df$rank[1],
      q_intercept = unname(stats::coef(model)[1]),
      q_slope = unname(stats::coef(model)[2]),
      q_r2 = summary(model)$r.squared
    )
  }))
  scaling <- merge(scaling, q_fit, by = "rank", all.x = TRUE, sort = FALSE)
  scaling$q_X_hat <- scaling$q_intercept + scaling$q_slope * scaling$U
  if (!is.null(sampling_rates)) {
    q_fit_rate <- do.call(rbind, lapply(split(scaling, scaling$rank), function(df) {
      model <- stats::lm(q_X ~ sampling_rate, data = df)
      data.frame(
        rank = df$rank[1],
        q_rate_intercept = unname(stats::coef(model)[1]),
        q_rate_slope = unname(stats::coef(model)[2]),
        q_rate_r2 = summary(model)$r.squared
      )
    }))
    scaling <- merge(scaling, q_fit_rate, by = "rank", all.x = TRUE, sort = FALSE)
    scaling$q_X_rate_hat <- scaling$q_rate_intercept + scaling$q_rate_slope * scaling$sampling_rate
  } else {
    q_fit_rate <- NULL
    scaling$q_X_rate_hat <- NA_real_
  }
  total_by <- c("U", "replicate")
  if (!is.null(sampling_rates)) total_by <- c("U", "sampling_rate", "replicate")
  q_total <- stats::aggregate(
    scaling["q_X"],
    by = scaling[total_by],
    FUN = sum
  )
  names(q_total)[names(q_total) == "q_X"] <- "q_total"
  total_model <- stats::lm(q_total ~ U, data = q_total)
  q_total$q_total_hat <- unname(stats::predict(total_model, newdata = q_total))
  q_total_fit <- data.frame(
    q_total_intercept = unname(stats::coef(total_model)[1]),
    q_total_slope = unname(stats::coef(total_model)[2]),
    q_total_r2 = summary(total_model)$r.squared
  )
  q_total_rate_fit <- NULL
  if (!is.null(sampling_rates)) {
    total_rate_model <- stats::lm(q_total ~ sampling_rate, data = q_total)
    q_total$q_total_rate_hat <- unname(stats::predict(total_rate_model, newdata = q_total))
    q_total_rate_fit <- data.frame(
      q_total_rate_intercept = unname(stats::coef(total_rate_model)[1]),
      q_total_rate_slope = unname(stats::coef(total_rate_model)[2]),
      q_total_rate_r2 = summary(total_rate_model)$r.squared
    )
  } else {
    q_total$q_total_rate_hat <- NA_real_
  }

  out <- list(
    scaling = scaling[order(scaling$U, scaling$replicate, scaling$rank), ],
    q_fit = q_fit[order(q_fit$rank), ],
    q_fit_rate = if (is.null(q_fit_rate)) NULL else q_fit_rate[order(q_fit_rate$rank), ],
    q_total = q_total[order(q_total$U, q_total$replicate), ],
    q_total_fit = q_total_fit,
    q_total_rate_fit = q_total_rate_fit,
    reference_fit = reference_fit,
    theta_Y = theta_Y,
    U_grid = U_grid,
    sampling_rates = sampling_rates,
    r = r_eff,
    n_replicates = n_replicates,
    empirical = isTRUE(empirical)
  )
  class(out) <- "scscale_umi_scaling"
  out
}

scscale_cell_scaling <- function(
  d2_X,
  p,
  n_grid,
  theta_Y,
  eps = 1e-12
) {
  if (inherits(d2_X, "scscale_fit")) {
    fit <- d2_X
    d2_X <- fit$spikes$d2_X
    p <- p %||% fit$p
  }
  if (missing(p) || is.null(p)) stop("p is required.", call. = FALSE)
  if (length(n_grid) < 1L) stop("n_grid must be non-empty.", call. = FALSE)

  rows <- lapply(n_grid, function(n) {
    c_X <- p / n
    q_X <- c_X * d2_X
    theta_X <- scscale_recoverability(q_X, c_X = c_X)
    data.frame(
      n = n,
      c_X = c_X,
      I_theory = scscale_information_from_theta(theta_X, theta_Y, eps = eps),
      rank = seq_along(q_X),
      d2_X = d2_X,
      q_X = q_X,
      theta_X = theta_X
    )
  })
  do.call(rbind, rows)
}
