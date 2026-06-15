# R-style fit objects for scScale models.

fit_cell_scaling <- function(
  counts,
  cell_grid,
  reps = 1,
  n_features = 300,
  transform = c("log1p", "pearson", "none"),
  min_cells = 10,
  r = 8,
  R = 3,
  p_sim = 100,
  seed = 1
) {
  transform <- match.arg(transform)
  counts <- counts_matrix(counts)
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

  cell_grid <- sort(unique(as.integer(cell_grid)))
  cell_grid <- cell_grid[is.finite(cell_grid) & cell_grid >= 5 & cell_grid <= ncol(counts)]
  if (length(cell_grid) < 3) stop("cell_grid must contain at least 3 feasible cell counts.", call. = FALSE)

  ref <- scaling_reference_fit(
    counts,
    n_features = n_features,
    transform = transform,
    min_cells = min_cells,
    r = r,
    R = R,
    p_sim = p_sim,
    seed = seed
  )

  rows <- list()
  idx <- 1L
  for (rep_id in seq_len(reps)) {
    for (n_cells in cell_grid) {
      set.seed(seed + 10000L * rep_id + n_cells)
      cells <- if (n_cells == ncol(counts)) colnames(counts) else sample(colnames(counts), n_cells)
      counts_sub <- counts[, cells, drop = FALSE]
      fit <- scaling_reference_fit(
        counts_sub,
        n_features = n_features,
        transform = transform,
        min_cells = min_cells,
        r = r,
        R = R,
        p_sim = p_sim,
        seed = seed
      )
      mi <- mi_theory(fit$spikes, ref$spikes, side = "cells", r = r)
      rows[[idx]] <- data.frame(
        rep = rep_id,
        n_cells = n_cells,
        mean_umi_per_cell = mean(colSums(as_dense_matrix(counts_sub))),
        noise_scale = fit$noise$scale,
        mp_upper = fit$spikes$mp_upper,
        n_spikes = fit$spikes$n_spikes,
        lambda1_over_mp_edge = max(fit$eigenvalues$eigenvalues) / fit$spikes$mp_upper,
        mi = mi$mi,
        mi_norm = mi$mi / r
      )
      idx <- idx + 1L
    }
  }

  observed <- do.call(rbind, rows)
  summary_df <- aggregate_scaling_observations(observed, x_col = "n_cells")
  nls_fit <- fit_saturating_nls(summary_df, x_col = "n_cells", y_col = "mean_mi_norm")
  new_scscale_fit(
    type = "cell",
    model = "saturating",
    x_col = "n_cells",
    y_col = "mean_mi_norm",
    fit = nls_fit,
    observed = observed,
    data = summary_df,
    reference = ref,
    settings = list(n_features = n_features, transform = transform, min_cells = min_cells, r = r, R = R, p_sim = p_sim, seed = seed),
    class = "scscale_cell_fit"
  )
}

fit_umi_scaling <- function(
  counts,
  umi_grid,
  n_cells = min(3000L, ncol(counts_matrix(counts))),
  reps = 1,
  n_features = 300,
  transform = c("log1p", "pearson", "none"),
  min_cells = 10,
  r = 8,
  R = 3,
  p_sim = 100,
  seed = 1
) {
  transform <- match.arg(transform)
  counts <- counts_matrix(counts)
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  n_cells <- min(as.integer(n_cells), ncol(counts))
  if (!is.finite(n_cells) || n_cells < 5) stop("n_cells must be at least 5.", call. = FALSE)

  umi_grid <- sort(unique(as.numeric(umi_grid)))
  umi_grid <- umi_grid[is.finite(umi_grid) & umi_grid > 0 & umi_grid <= 1]
  if (length(umi_grid) < 3) stop("umi_grid must contain at least 3 fractions in (0, 1].", call. = FALSE)

  set.seed(seed)
  cells_base <- if (n_cells == ncol(counts)) colnames(counts) else sample(colnames(counts), n_cells)
  counts_base <- counts[, cells_base, drop = FALSE]
  ref <- scaling_reference_fit(
    counts_base,
    n_features = n_features,
    transform = transform,
    min_cells = min_cells,
    r = r,
    R = R,
    p_sim = p_sim,
    seed = seed
  )

  rows <- list()
  idx <- 1L
  for (rep_id in seq_len(reps)) {
    for (frac in umi_grid) {
      set.seed(seed + 10000L * rep_id + as.integer(round(frac * 10000)))
      counts_sub <- thin_counts_fraction(counts_base, frac)
      fit <- scaling_reference_fit(
        counts_sub,
        n_features = n_features,
        transform = transform,
        min_cells = min_cells,
        r = r,
        R = R,
        p_sim = p_sim,
        seed = seed
      )
      mi <- mi_theory(fit$spikes, ref$spikes, side = "cells", r = r)
      rows[[idx]] <- data.frame(
        rep = rep_id,
        umi_fraction = frac,
        mean_umi_per_cell = mean(colSums(as_dense_matrix(counts_sub))),
        noise_scale = fit$noise$scale,
        mp_upper = fit$spikes$mp_upper,
        n_spikes = fit$spikes$n_spikes,
        lambda1_over_mp_edge = max(fit$eigenvalues$eigenvalues) / fit$spikes$mp_upper,
        mi = mi$mi,
        mi_norm = mi$mi / r
      )
      idx <- idx + 1L
    }
  }

  observed <- do.call(rbind, rows)
  summary_df <- aggregate_scaling_observations(observed, x_col = "mean_umi_per_cell")
  nls_fit <- fit_saturating_nls(summary_df, x_col = "mean_umi_per_cell", y_col = "mean_mi_norm")
  new_scscale_fit(
    type = "umi",
    model = "saturating",
    x_col = "mean_umi_per_cell",
    y_col = "mean_mi_norm",
    fit = nls_fit,
    observed = observed,
    data = summary_df,
    reference = ref,
    settings = list(umi_grid = umi_grid, n_cells = n_cells, n_features = n_features, transform = transform, min_cells = min_cells, r = r, R = R, p_sim = p_sim, seed = seed),
    class = "scscale_umi_fit"
  )
}

fit_batch_scaling <- function(
  summary_df,
  law = c("batch_number", "cells_per_batch"),
  target_col = "mean_I_bio_norm",
  min_points = 5
) {
  law <- match.arg(law)
  fit <- fit_batch_effect_scaling(
    summary_df = summary_df,
    law = law,
    target_col = target_col,
    min_points = min_points
  )
  x_col <- if (law == "batch_number") "m_batch" else "cells_per_batch"
  class(fit) <- c("scscale_batch_fit", "scscale_fit", class(fit))
  fit$type <- "batch"
  fit$model <- law
  fit$x_col <- x_col
  fit$y_col <- "I_fit"
  fit$observed <- fit$data
  fit
}

predict.scscale_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$pred$I_pred)
  dat <- as.data.frame(newdata)
  if (!object$x_col %in% colnames(dat)) {
    stop("newdata must contain column '", object$x_col, "'.", call. = FALSE)
  }
  x <- dat[[object$x_col]]
  co <- object$coef

  if (identical(object$model, "saturating")) {
    return(unname(co[["I_inf"]] * (1 - exp(-co[["k"]] * x))))
  }
  if (identical(object$model, "batch_number")) {
    return(unname(co[["I_inf"]] - co[["C"]] * log(1 - co[["A"]] / x)))
  }
  if (identical(object$model, "cells_per_batch")) {
    return(unname(co[["I_inf"]] - co[["C"]] * log(1 + co[["A"]] / x)))
  }
  stop("Unknown scScale model: ", object$model, call. = FALSE)
}

coef.scscale_fit <- function(object, ...) {
  object$coef
}

summary.scscale_fit <- function(object, ...) {
  out <- data.frame(
    type = object$type,
    model = object$model,
    x_col = object$x_col,
    y_col = object$y_col,
    n_points = nrow(object$data),
    ok = isTRUE(object$ok %||% TRUE),
    message = object$message %||% "ok"
  )
  cbind(out, as.data.frame(as.list(object$coef)), object$stats)
}

plot.scscale_fit <- function(x, xlab = x$x_col, ylab = x$y_col, ...) {
  y <- if (x$y_col %in% colnames(x$data)) x$data[[x$y_col]] else x$data$mean_mi_norm
  graphics::plot(
    x$data[[x$x_col]],
    y,
    pch = 19,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  pred <- scscale_prediction_grid(x)
  graphics::lines(pred[[x$x_col]], pred$I_pred, col = "#1f77b4", lwd = 2)
  invisible(x)
}

new_scscale_fit <- function(type, model, x_col, y_col, fit, observed, data, reference, settings, class) {
  pred <- scscale_prediction_grid_from_fit(fit, model = model, x_col = x_col, data = data)
  data$I_pred <- scscale_model_predict(fit, model = model, x_col = x_col, newdata = data)
  data$resid <- data[[y_col]] - data$I_pred
  structure(
    list(
      ok = TRUE,
      type = type,
      model = model,
      x_col = x_col,
      y_col = y_col,
      fit = fit,
      coef = scscale_model_coef(fit),
      stats = batch_fit_stats(data[[y_col]], data$I_pred),
      observed = observed,
      data = data,
      pred = pred,
      reference = reference,
      settings = settings,
      message = "ok"
    ),
    class = c(class, "scscale_fit")
  )
}

scaling_reference_fit <- function(counts, n_features, transform, min_cells, r, R, p_sim, seed) {
  ev <- find_eigenvalues(
    counts,
    n_features = n_features,
    transform = transform,
    min_cells = min_cells
  )
  noise <- fit_noise(ev, R = R, p_sim = p_sim, seed = seed)
  spikes <- fit_spikes(ev, noise = noise, r = r)
  list(eigenvalues = ev, noise = noise, spikes = spikes)
}

thin_counts_fraction <- function(counts, fraction) {
  X <- as_dense_matrix(counts)
  matrix(
    stats::rbinom(length(X), size = as.integer(X), prob = fraction),
    nrow = nrow(X),
    dimnames = dimnames(X)
  )
}

aggregate_scaling_observations <- function(observed, x_col) {
  split_key <- observed[[x_col]]
  rows <- lapply(split(observed, split_key), function(x) {
    out <- x[1, x_col, drop = FALSE]
    out$mean_mi <- mean(x$mi, na.rm = TRUE)
    out$sd_mi <- stats::sd(x$mi, na.rm = TRUE)
    out$se_mi <- out$sd_mi / sqrt(nrow(x))
    out$mean_mi_norm <- mean(x$mi_norm, na.rm = TRUE)
    out$sd_mi_norm <- stats::sd(x$mi_norm, na.rm = TRUE)
    out$se_mi_norm <- out$sd_mi_norm / sqrt(nrow(x))
    out$mean_lambda1_over_mp_edge <- mean(x$lambda1_over_mp_edge, na.rm = TRUE)
    out$mean_n_spikes <- mean(x$n_spikes, na.rm = TRUE)
    out$n_rep_observed <- nrow(x)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out[[x_col]]), , drop = FALSE]
}

fit_saturating_nls <- function(dat, x_col, y_col) {
  fit_dat <- dat
  fit_dat$x <- fit_dat[[x_col]]
  fit_dat$y <- fit_dat[[y_col]]
  y_max <- max(fit_dat$y, na.rm = TRUE)
  x_med <- stats::median(fit_dat$x)
  starts <- list(
    list(I_inf = y_max, k = 1 / x_med),
    list(I_inf = y_max * 1.05 + 1e-8, k = 0.1 / x_med),
    list(I_inf = y_max * 1.25 + 1e-8, k = 2 / x_med),
    list(I_inf = mean(fit_dat$y, na.rm = TRUE) + 1e-8, k = 1 / max(fit_dat$x))
  )
  best <- NULL
  best_rmse <- Inf
  for (start in starts) {
    fit <- tryCatch(
      stats::nls(
        y ~ I_inf * (1 - exp(-k * x)),
        data = fit_dat,
        start = start,
        algorithm = "port",
        lower = c(I_inf = 0, k = 0),
        control = stats::nls.control(maxiter = 1000, warnOnly = TRUE)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      pred <- as.numeric(stats::predict(fit, newdata = fit_dat))
      rmse <- sqrt(mean((fit_dat$y - pred)^2))
      if (is.finite(rmse) && rmse < best_rmse) {
        best <- fit
        best_rmse <- rmse
      }
    }
  }
  if (!is.null(best)) return(best)

  # Flat or nearly flat curves can make the nonlinear gradient singular.
  # Use a stable near-flat saturating fit so predict()/plot() still work.
  structure(
    list(coef = c(I_inf = mean(fit_dat$y, na.rm = TRUE), k = 0)),
    class = "scscale_manual_saturating_fit"
  )
}

scscale_prediction_grid <- function(object) {
  scscale_prediction_grid_from_fit(object$fit, object$model, object$x_col, object$data)
}

scscale_prediction_grid_from_fit <- function(fit, model, x_col, data) {
  x <- seq(min(data[[x_col]]), max(data[[x_col]]), length.out = 200)
  pred <- data.frame(x = x)
  names(pred) <- x_col
  co <- scscale_model_coef(fit)
  if (identical(model, "saturating")) {
    pred$I_pred <- co[["I_inf"]] * (1 - exp(-co[["k"]] * pred[[x_col]]))
  } else if (identical(model, "batch_number")) {
    pred$I_pred <- co[["I_inf"]] - co[["C"]] * log(1 - co[["A"]] / pred[[x_col]])
  } else if (identical(model, "cells_per_batch")) {
    pred$I_pred <- co[["I_inf"]] - co[["C"]] * log(1 + co[["A"]] / pred[[x_col]])
  }
  pred
}

scscale_model_coef <- function(fit) {
  if (inherits(fit, "scscale_manual_saturating_fit")) return(fit$coef)
  stats::coef(fit)
}

scscale_model_predict <- function(fit, model, x_col, newdata) {
  if (!inherits(fit, "scscale_manual_saturating_fit")) {
    return(as.numeric(stats::predict(fit, newdata = newdata)))
  }
  co <- scscale_model_coef(fit)
  x <- newdata[[x_col]]
  if (identical(model, "saturating")) {
    return(unname(co[["I_inf"]] * (1 - exp(-co[["k"]] * x))))
  }
  stop("Manual fit is only implemented for saturating models.", call. = FALSE)
}
