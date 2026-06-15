# Batch-effect scaling utilities adapted from the BMMC CITE-seq validation
# workflow in ChesleaZ/Scaling-law-summer-version-.

batch_subspace_mi <- function(Z1, Z2, eps = 1e-12) {
  common <- intersect(rownames(Z1), rownames(Z2))
  if (length(common) < 5) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  Z1 <- as.matrix(Z1[common, , drop = FALSE])
  Z2 <- as.matrix(Z2[common, , drop = FALSE])
  Q1 <- orthonormalize_subspace(Z1)
  Q2 <- orthonormalize_subspace(Z2)
  if (is.null(Q1) || is.null(Q2)) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  r_eff <- min(ncol(Q1), ncol(Q2))
  if (r_eff < 1) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  gamma <- svd(crossprod(Q1, Q2), nu = 0, nv = 0)$d^2
  gamma <- gamma[seq_len(min(r_eff, length(gamma)))]
  gamma <- pmin(pmax(gamma, 0), 1 - eps)
  mi <- overlap_mi_from_gamma(gamma)

  list(mi = mi, mi_norm = mi / length(gamma), gamma = gamma, r_eff = length(gamma))
}

batch_effect_mi <- function(
  x,
  y,
  r_x = 20,
  r_y = 20,
  transform_x = c("log1p", "pearson", "none"),
  transform_y = c("log1p", "clr", "none"),
  n_features = NULL,
  min_cells = 10,
  clip = 10
) {
  transform_x <- match.arg(transform_x)
  transform_y <- match.arg(transform_y)
  x <- counts_matrix(x)
  y <- counts_matrix(y)
  if (is.null(colnames(x)) && is.null(colnames(y)) && ncol(x) == ncol(y)) {
    cell_names <- paste0("cell_", seq_len(ncol(x)))
    colnames(x) <- cell_names
    colnames(y) <- cell_names
  } else {
    if (is.null(colnames(x))) colnames(x) <- paste0("x_cell_", seq_len(ncol(x)))
    if (is.null(colnames(y))) colnames(y) <- paste0("y_cell_", seq_len(ncol(y)))
  }

  common <- intersect(colnames(x), colnames(y))
  if (length(common) < 5) {
    stop("x and y must share at least 5 cells by column name.", call. = FALSE)
  }
  x <- x[, common, drop = FALSE]
  y <- y[, common, drop = FALSE]

  if (!is.null(n_features)) {
    features <- select_hvgs(x, n_features = n_features, min_cells = min_cells)
    x <- x[features, , drop = FALSE]
  }

  X <- transform_batch_matrix(x, transform_x, clip = clip)
  Y <- transform_batch_matrix(y, transform_y, clip = clip)
  Zx <- right_singular_vectors(X, r = r_x)
  Zy <- right_singular_vectors(Y, r = r_y)
  if (is.null(Zx) || is.null(Zy)) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0L))
  }

  batch_subspace_mi(Zx, Zy)
}

sample_batch_cells <- function(meta, batch_col, m_batch, cells_per_batch, seed = NULL, batches = NULL) {
  if (!batch_col %in% colnames(meta)) {
    stop("batch_col is not a column in meta.", call. = FALSE)
  }
  if (is.null(rownames(meta))) {
    stop("meta must have cell names in rownames.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  batch_values <- as.character(meta[[batch_col]])
  sizes <- table(batch_values)
  eligible <- names(sizes)[sizes >= cells_per_batch]
  if (is.null(batches)) {
    if (length(eligible) < m_batch) return(NULL)
    batches <- sample(eligible, m_batch, replace = FALSE)
  } else {
    batches <- as.character(batches)
    if (length(batches) != m_batch || any(!batches %in% eligible)) return(NULL)
  }

  cells <- unlist(lapply(batches, function(batch) {
    sample(rownames(meta)[batch_values == batch], cells_per_batch, replace = FALSE)
  }), use.names = FALSE)

  list(cells = cells, batches = batches)
}

fit_batch_effect_scaling <- function(
  summary_df,
  law = c("batch_number", "cells_per_batch"),
  target_col = "mean_I_bio_norm",
  min_points = 5
) {
  law <- match.arg(law)
  dat <- as.data.frame(summary_df)
  x_col <- if (law == "batch_number") "m_batch" else "cells_per_batch"
  if (!all(c(x_col, target_col) %in% colnames(dat))) {
    stop("summary_df must contain ", x_col, " and ", target_col, ".", call. = FALSE)
  }

  dat$I_fit <- dat[[target_col]]
  dat <- dat[is.finite(dat$I_fit) & is.finite(dat[[x_col]]) & dat[[x_col]] > 0, , drop = FALSE]
  dat <- dat[order(dat[[x_col]]), , drop = FALSE]

  if (nrow(dat) < min_points) {
    return(batch_scaling_failure(dat, paste0("Too few points for fit: ", nrow(dat), " < ", min_points)))
  }

  if (law == "batch_number") {
    fit <- fit_batch_number_nls(dat, x_col = x_col)
    pred_x <- seq(min(dat[[x_col]]), max(dat[[x_col]]), length.out = 200)
    pred <- data.frame(m_batch = pred_x)
    if (!is.null(fit$fit)) {
      co <- stats::coef(fit$fit)
      pred$I_pred <- co[["I_inf"]] - co[["C"]] * log(1 - co[["A"]] / pred$m_batch)
    }
  } else {
    fit <- fit_cells_per_batch_nls(dat, x_col = x_col)
    pred_x <- seq(min(dat[[x_col]]), max(dat[[x_col]]), length.out = 200)
    pred <- data.frame(cells_per_batch = pred_x)
    if (!is.null(fit$fit)) {
      co <- stats::coef(fit$fit)
      pred$I_pred <- co[["I_inf"]] - co[["C"]] * log(1 + co[["A"]] / pred$cells_per_batch)
    }
  }

  if (is.null(fit$fit)) {
    return(batch_scaling_failure(dat, fit$error))
  }

  dat$I_pred <- as.numeric(stats::predict(fit$fit, newdata = dat))
  dat$resid <- dat$I_fit - dat$I_pred
  structure(
    list(
      ok = TRUE,
      law = law,
      target_col = target_col,
      fit = fit$fit,
      data = dat,
      pred = pred,
      coef = stats::coef(fit$fit),
      stats = batch_fit_stats(dat$I_fit, dat$I_pred),
      message = "ok"
    ),
    class = "scsm_batch_scaling"
  )
}

summarize_batch_effect_results <- function(results) {
  dat <- as.data.frame(results)
  if (!"status" %in% colnames(dat)) dat$status <- "ok"
  dat <- dat[dat$status == "ok", , drop = FALSE]
  if (nrow(dat) == 0) return(data.frame())

  keys <- intersect(c("setting_id", "experiment", "m_batch", "cells_per_batch", "n_cells"), colnames(dat))
  if (length(keys) == 0) {
    stop("results must contain at least one grouping column.", call. = FALSE)
  }
  split_key <- interaction(dat[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(dat, split_key), function(x) {
    out <- x[1, keys, drop = FALSE]
    out$mean_I_bio <- mean(x$I_bio, na.rm = TRUE)
    out$sd_I_bio <- stats::sd(x$I_bio, na.rm = TRUE)
    out$se_I_bio <- out$sd_I_bio / sqrt(nrow(x))
    out$mean_I_bio_norm <- mean(x$I_bio_norm, na.rm = TRUE)
    out$sd_I_bio_norm <- stats::sd(x$I_bio_norm, na.rm = TRUE)
    out$se_I_bio_norm <- out$sd_I_bio_norm / sqrt(nrow(x))
    out$n_rep_observed <- nrow(x)
    out
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  order_cols <- intersect(c("setting_id", "experiment", "m_batch", "cells_per_batch"), colnames(out))
  if (length(order_cols) > 0) {
    out <- out[do.call(order, out[order_cols]), , drop = FALSE]
  }
  out
}

transform_batch_matrix <- function(x, transform, clip = 10) {
  if (transform == "pearson") return(pearson_residuals(x, clip = clip))
  X <- as_dense_matrix(x)
  if (transform == "log1p") {
    X <- log1p(X)
  } else if (transform == "clr") {
    denom <- exp(colMeans(log1p(X)))
    denom[!is.finite(denom) | denom <= 0] <- 1
    X <- log1p(t(t(X) / denom))
  }
  X <- t(scale(t(X), center = TRUE, scale = TRUE))
  X[!is.finite(X)] <- 0
  rownames(X) <- rownames(x)
  colnames(X) <- colnames(x)
  X
}

right_singular_vectors <- function(x, r = 20) {
  X <- as.matrix(x)
  r_eff <- min(r, nrow(X) - 1L, ncol(X) - 1L)
  if (r_eff < 1) return(NULL)

  if (requireNamespace("irlba", quietly = TRUE) && r_eff < min(dim(X)) - 1L) {
    sv <- irlba::irlba(X, nv = r_eff, nu = 0, maxit = 2000, work = max(3L * r_eff, 30L))
    Z <- sv$v[, seq_len(r_eff), drop = FALSE]
  } else {
    sv <- svd(X, nu = 0, nv = r_eff)
    Z <- sv$v[, seq_len(r_eff), drop = FALSE]
  }
  rownames(Z) <- colnames(X)
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  Z
}

orthonormalize_subspace <- function(Z) {
  Z <- scale(as.matrix(Z), center = TRUE, scale = FALSE)
  Z[!is.finite(Z)] <- 0
  keep <- apply(Z, 2, stats::sd) > 1e-12
  Z <- Z[, keep, drop = FALSE]
  if (ncol(Z) < 1) return(NULL)

  qrZ <- qr(Z)
  if (qrZ$rank < 1) return(NULL)
  Q <- qr.Q(qrZ)[, seq_len(qrZ$rank), drop = FALSE]
  rownames(Q) <- rownames(Z)
  Q
}

batch_fit_stats <- function(y, pred) {
  resid <- y - pred
  denom <- sum((y - mean(y))^2)
  data.frame(
    R2 = if (denom > 0) 1 - sum(resid^2) / denom else NA_real_,
    RMSE = sqrt(mean(resid^2)),
    MAE = mean(abs(resid))
  )
}

try_batch_nls <- function(expr, dat, starts, lower, upper) {
  best <- NULL
  best_rmse <- Inf
  errors <- character(0)
  for (start in starts) {
    fit <- tryCatch(
      stats::nls(
        expr,
        data = dat,
        start = start,
        algorithm = "port",
        lower = lower,
        upper = upper,
        control = stats::nls.control(maxiter = 1000, warnOnly = TRUE)
      ),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      }
    )
    if (!is.null(fit)) {
      pred <- as.numeric(stats::predict(fit, newdata = dat))
      rmse <- sqrt(mean((dat$I_fit - pred)^2))
      if (is.finite(rmse) && rmse < best_rmse) {
        best <- fit
        best_rmse <- rmse
      }
    }
  }
  if (is.null(best)) return(list(fit = NULL, error = paste(unique(errors), collapse = " | ")))
  list(fit = best, error = NA_character_)
}

fit_batch_number_nls <- function(dat, x_col) {
  names(dat)[names(dat) == x_col] <- "m_batch"
  min_m <- min(dat$m_batch)
  max_abs_a <- 0.95 * min_m
  y_rng <- diff(range(dat$I_fit))
  c_start <- max(abs(y_rng), stats::sd(dat$I_fit), 1e-3)
  starts <- list(
    list(I_inf = mean(dat$I_fit), C = c_start, A = 0.25 * min_m),
    list(I_inf = min(dat$I_fit), C = c_start, A = -0.25 * min_m),
    list(I_inf = max(dat$I_fit), C = 0.5 * c_start, A = 0.5 * min_m),
    list(I_inf = mean(dat$I_fit), C = 2 * c_start, A = -0.5 * min_m)
  )
  try_batch_nls(
    I_fit ~ I_inf - C * log(1 - A / m_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = -max_abs_a),
    upper = c(I_inf = Inf, C = Inf, A = max_abs_a)
  )
}

fit_cells_per_batch_nls <- function(dat, x_col) {
  names(dat)[names(dat) == x_col] <- "cells_per_batch"
  y_rng <- diff(range(dat$I_fit))
  c_start <- max(abs(y_rng), stats::sd(dat$I_fit), 1e-3)
  s_med <- stats::median(dat$cells_per_batch)
  s_max <- max(dat$cells_per_batch)
  starts <- list(
    list(I_inf = max(dat$I_fit), C = c_start, A = s_med),
    list(I_inf = mean(dat$I_fit), C = 0.5 * c_start, A = 0.5 * s_med),
    list(I_inf = max(dat$I_fit) + c_start, C = 2 * c_start, A = 2 * s_med),
    list(I_inf = min(dat$I_fit), C = c_start, A = 0.1 * s_med),
    list(I_inf = max(dat$I_fit), C = 0.25 * c_start, A = 5 * s_med),
    list(I_inf = max(dat$I_fit) + 0.5 * c_start, C = 4 * c_start, A = 0.25 * s_med)
  )
  try_batch_nls(
    I_fit ~ I_inf - C * log(1 + A / cells_per_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = 0),
    upper = c(I_inf = Inf, C = Inf, A = 20 * s_max)
  )
}

batch_scaling_failure <- function(dat, message) {
  structure(
    list(
      ok = FALSE,
      data = dat,
      pred = data.frame(),
      message = message
    ),
    class = "scsm_batch_scaling"
  )
}
