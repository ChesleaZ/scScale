# Direct empirical scaling-law fits.

scscale_empirical_cell_curve <- function(n, I_infinity, C_n, A_n) {
  I_infinity - C_n * log1p(A_n / n)
}

scscale_empirical_umi_curve <- function(U, I_infinity, C_U, A_U, U_offset) {
  I_infinity - C_U * log1p(A_U / (U + U_offset))
}

scscale_empirical_joint_curve <- function(n, U, I_infinity, C_n, A_n, C_U, A_U, U_offset) {
  I_infinity -
    C_n * log1p(A_n / n) -
    C_U * log1p(A_U / (U + U_offset))
}

scscale_empirical_gap_cell <- function(n, C_n, A_n) {
  C_n * log1p(A_n / n)
}

scscale_empirical_gap_umi <- function(U, C_U, A_U, U_offset) {
  C_U * log1p(A_U / (U + U_offset))
}

scscale_scaling_fit_metrics <- function(y, fitted, k) {
  residual <- y - fitted
  rss <- sum(residual^2)
  n <- length(y)
  data.frame(
    n = n,
    df = max(n - k, 0L),
    rss = rss,
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    aic = if (n > 0L && rss > 0) n * log(rss / n) + 2 * k else NA_real_
  )
}

scscale_clean_scaling_data <- function(data, x_cols, y_col) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  needed <- c(x_cols, y_col)
  missing <- setdiff(needed, names(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data[, needed, drop = FALSE]
  for (nm in needed) out[[nm]] <- as.numeric(out[[nm]])
  keep <- stats::complete.cases(out)
  for (nm in x_cols) keep <- keep & is.finite(out[[nm]]) & out[[nm]] > 0
  keep <- keep & is.finite(out[[y_col]])
  out <- out[keep, , drop = FALSE]
  if (nrow(out) < 3L) stop("at least three complete positive observations are required.", call. = FALSE)
  out
}

scscale_positive_start <- function(x, fallback = 1) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x)) stats::median(x) else fallback
}

scscale_fit_empirical_cell_law <- function(
  data,
  n_col = "n",
  mi_col = "I_empirical",
  start = NULL,
  control = list(maxit = 10000)
) {
  df <- scscale_clean_scaling_data(data, n_col, mi_col)
  n <- df[[n_col]]
  y <- df[[mi_col]]
  ceiling_min <- max(y)
  n_scale <- scscale_positive_start(n)
  y_span <- max(diff(range(y)), stats::sd(y), 1e-3, na.rm = TRUE)
  start <- utils::modifyList(
    list(I_extra = max(y_span, 1e-3), C_n = max(y_span, 1e-3), A_n = n_scale),
    start %||% list()
  )
  par0 <- log(c(start$I_extra, start$C_n, start$A_n))

  unpack <- function(par) {
    vals <- exp(par)
    list(I_infinity = ceiling_min + vals[[1]], C_n = vals[[2]], A_n = vals[[3]])
  }
  objective <- function(par) {
    p <- unpack(par)
    fitted <- scscale_empirical_cell_curve(n, p$I_infinity, p$C_n, p$A_n)
    sum((y - fitted)^2)
  }
  opt <- stats::optim(par0, objective, method = "Nelder-Mead", control = control)
  p <- unpack(opt$par)
  fitted <- scscale_empirical_cell_curve(n, p$I_infinity, p$C_n, p$A_n)
  coef <- c(I_infinity = p$I_infinity, C_n = p$C_n, A_n = p$A_n)
  out <- list(
    model = "cell",
    formula = "I(n) = I_infinity - C_n * log(1 + A_n / n)",
    coefficients = coef,
    fixed = numeric(0),
    data = df,
    fitted = fitted,
    residuals = y - fitted,
    metrics = scscale_scaling_fit_metrics(y, fitted, length(coef)),
    convergence = opt$convergence,
    objective = opt$value,
    columns = list(n = n_col, mi = mi_col),
    ceiling_min = ceiling_min
  )
  class(out) <- "scscale_empirical_scaling_fit"
  out
}

scscale_fit_empirical_umi_law <- function(
  data,
  U_col = "U",
  mi_col = "I_empirical",
  start = NULL,
  control = list(maxit = 10000)
) {
  df <- scscale_clean_scaling_data(data, U_col, mi_col)
  U <- df[[U_col]]
  y <- df[[mi_col]]
  ceiling_min <- max(y)
  U_scale <- scscale_positive_start(U)
  y_span <- max(diff(range(y)), stats::sd(y), 1e-3, na.rm = TRUE)
  start <- utils::modifyList(
    list(I_extra = max(y_span, 1e-3), C_U = max(y_span, 1e-3), A_U = U_scale, U_offset = U_scale / 10),
    start %||% list()
  )
  par0 <- log(c(start$I_extra, start$C_U, start$A_U, start$U_offset))

  unpack <- function(par) {
    vals <- exp(par)
    list(I_infinity = ceiling_min + vals[[1]], C_U = vals[[2]], A_U = vals[[3]], U_offset = vals[[4]])
  }
  objective <- function(par) {
    p <- unpack(par)
    fitted <- scscale_empirical_umi_curve(U, p$I_infinity, p$C_U, p$A_U, p$U_offset)
    sum((y - fitted)^2)
  }
  opt <- stats::optim(par0, objective, method = "Nelder-Mead", control = control)
  p <- unpack(opt$par)
  fitted <- scscale_empirical_umi_curve(U, p$I_infinity, p$C_U, p$A_U, p$U_offset)
  coef <- c(I_infinity = p$I_infinity, C_U = p$C_U, A_U = p$A_U, U_offset = p$U_offset)
  out <- list(
    model = "umi",
    formula = "I(U) = I_infinity - C_U * log(1 + A_U / (U + U_offset))",
    coefficients = coef,
    fixed = numeric(0),
    data = df,
    fitted = fitted,
    residuals = y - fitted,
    metrics = scscale_scaling_fit_metrics(y, fitted, length(coef)),
    convergence = opt$convergence,
    objective = opt$value,
    columns = list(U = U_col, mi = mi_col),
    ceiling_min = ceiling_min
  )
  class(out) <- "scscale_empirical_scaling_fit"
  out
}

scscale_fit_empirical_joint_law <- function(
  data,
  n_col = "n",
  U_col = "U",
  mi_col = "I_empirical",
  cell_fit = NULL,
  umi_fit = NULL,
  mode = c("hybrid", "transfer", "free"),
  start = NULL,
  control = list(maxit = 20000)
) {
  mode <- match.arg(mode)
  df <- scscale_clean_scaling_data(data, c(n_col, U_col), mi_col)
  n <- df[[n_col]]
  U <- df[[U_col]]
  y <- df[[mi_col]]
  ceiling_min <- max(y)
  y_span <- max(diff(range(y)), stats::sd(y), 1e-3, na.rm = TRUE)

  cell_coef <- if (!is.null(cell_fit)) cell_fit$coefficients else NULL
  umi_coef <- if (!is.null(umi_fit)) umi_fit$coefficients else NULL
  start <- utils::modifyList(
    list(
      I_extra = max(y_span, 1e-3),
      C_n = cell_coef[["C_n"]] %||% max(y_span / 2, 1e-3),
      A_n = cell_coef[["A_n"]] %||% scscale_positive_start(n),
      C_U = umi_coef[["C_U"]] %||% max(y_span / 2, 1e-3),
      A_U = umi_coef[["A_U"]] %||% scscale_positive_start(U),
      U_offset = umi_coef[["U_offset"]] %||% (scscale_positive_start(U) / 10)
    ),
    start %||% list()
  )

  if (mode == "transfer") {
    needed <- c("C_n", "A_n", "C_U", "A_U", "U_offset")
    missing <- needed[!is.finite(unlist(start[needed])) | unlist(start[needed]) <= 0]
    if (length(missing)) {
      stop("transfer mode requires positive finite ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    gaps <- scscale_empirical_gap_cell(n, start$C_n, start$A_n) +
      scscale_empirical_gap_umi(U, start$C_U, start$A_U, start$U_offset)
    I_infinity <- max(ceiling_min, mean(y + gaps))
    fitted <- I_infinity - gaps
    coef <- c(I_infinity = I_infinity)
    fixed <- c(C_n = start$C_n, A_n = start$A_n, C_U = start$C_U, A_U = start$A_U, U_offset = start$U_offset)
    out <- list(
      model = "joint",
      mode = mode,
      formula = "I(n,U) = I_infinity - C_n * log(1 + A_n / n) - C_U * log(1 + A_U / (U + U_offset))",
      coefficients = coef,
      fixed = fixed,
      data = df,
      fitted = fitted,
      residuals = y - fitted,
      metrics = scscale_scaling_fit_metrics(y, fitted, length(coef)),
      convergence = 0L,
      objective = sum((y - fitted)^2),
      columns = list(n = n_col, U = U_col, mi = mi_col),
      ceiling_min = ceiling_min
    )
    class(out) <- "scscale_empirical_scaling_fit"
    return(out)
  }

  if (mode == "hybrid") {
    fixed <- c(A_n = start$A_n, A_U = start$A_U, U_offset = start$U_offset)
    par0 <- log(c(start$I_extra, start$C_n, start$C_U))
    unpack <- function(par) {
      vals <- exp(par)
      c(I_infinity = ceiling_min + vals[[1]], C_n = vals[[2]], C_U = vals[[3]], fixed)
    }
  } else {
    fixed <- numeric(0)
    par0 <- log(c(start$I_extra, start$C_n, start$A_n, start$C_U, start$A_U, start$U_offset))
    unpack <- function(par) {
      vals <- exp(par)
      c(
        I_infinity = ceiling_min + vals[[1]],
        C_n = vals[[2]],
        A_n = vals[[3]],
        C_U = vals[[4]],
        A_U = vals[[5]],
        U_offset = vals[[6]]
      )
    }
  }
  objective <- function(par) {
    p <- unpack(par)
    fitted <- scscale_empirical_joint_curve(n, U, p[["I_infinity"]], p[["C_n"]], p[["A_n"]], p[["C_U"]], p[["A_U"]], p[["U_offset"]])
    sum((y - fitted)^2)
  }
  opt <- stats::optim(par0, objective, method = "Nelder-Mead", control = control)
  coef <- unpack(opt$par)
  fitted <- scscale_empirical_joint_curve(n, U, coef[["I_infinity"]], coef[["C_n"]], coef[["A_n"]], coef[["C_U"]], coef[["A_U"]], coef[["U_offset"]])
  out <- list(
    model = "joint",
    mode = mode,
    formula = "I(n,U) = I_infinity - C_n * log(1 + A_n / n) - C_U * log(1 + A_U / (U + U_offset))",
    coefficients = coef[setdiff(names(coef), names(fixed))],
    fixed = fixed,
    data = df,
    fitted = fitted,
    residuals = y - fitted,
    metrics = scscale_scaling_fit_metrics(y, fitted, length(coef) - length(fixed)),
    convergence = opt$convergence,
    objective = opt$value,
    columns = list(n = n_col, U = U_col, mi = mi_col),
    ceiling_min = ceiling_min
  )
  class(out) <- "scscale_empirical_scaling_fit"
  out
}

scscale_empirical_global_ceiling <- function(cell_fit, umi_fit, n0, U0) {
  if (!inherits(cell_fit, "scscale_empirical_scaling_fit") || cell_fit$model != "cell") {
    stop("cell_fit must be a cell empirical scaling fit.", call. = FALSE)
  }
  if (!inherits(umi_fit, "scscale_empirical_scaling_fit") || umi_fit$model != "umi") {
    stop("umi_fit must be a UMI empirical scaling fit.", call. = FALSE)
  }
  cc <- cell_fit$coefficients
  uc <- umi_fit$coefficients
  from_cell <- cc[["I_infinity"]] + scscale_empirical_gap_umi(U0, uc[["C_U"]], uc[["A_U"]], uc[["U_offset"]])
  from_umi <- uc[["I_infinity"]] + scscale_empirical_gap_cell(n0, cc[["C_n"]], cc[["A_n"]])
  c(from_cell_marginal = from_cell, from_umi_marginal = from_umi, average = mean(c(from_cell, from_umi)))
}

scscale_nearest_grid_value <- function(values, target, label) {
  values <- sort(unique(as.numeric(values)))
  values <- values[is.finite(values)]
  if (!length(values)) stop(label, " grid is empty.", call. = FALSE)
  if (is.null(target)) return(max(values))
  target <- as.numeric(target)
  if (!is.finite(target)) stop(label, " reference must be finite.", call. = FALSE)
  values[which.min(abs(values - target))]
}

scscale_fit_empirical_laws_from_grid <- function(
  grid,
  n_ref = NULL,
  U_ref = NULL,
  joint_mode = c("hybrid", "transfer", "free"),
  cell_start = NULL,
  umi_start = NULL,
  joint_start = NULL,
  cell_control = list(maxit = 10000),
  umi_control = list(maxit = 10000),
  joint_control = list(maxit = 20000)
) {
  joint_mode <- match.arg(joint_mode)
  if (!is.data.frame(grid)) grid <- as.data.frame(grid)
  needed <- c("n", "cell_sampling_rate", "U", "sampling_rate", "I_empirical")
  missing <- setdiff(needed, names(grid))
  if (length(missing)) {
    stop("grid is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  n_ref <- scscale_nearest_grid_value(grid$n, n_ref, "n")
  U_ref <- scscale_nearest_grid_value(grid$U, U_ref, "U")
  cell_data <- grid[grid$U == U_ref, , drop = FALSE]
  umi_data <- grid[grid$n == n_ref, , drop = FALSE]
  if (nrow(cell_data) < 3L) {
    stop("cell marginal fit needs at least three grid rows at U_ref = ", U_ref, ".", call. = FALSE)
  }
  if (nrow(umi_data) < 3L) {
    stop("UMI marginal fit needs at least three grid rows at n_ref = ", n_ref, ".", call. = FALSE)
  }

  cell_fit <- scscale_fit_empirical_cell_law(
    cell_data,
    n_col = "n",
    mi_col = "I_empirical",
    start = cell_start,
    control = cell_control
  )
  umi_fit <- scscale_fit_empirical_umi_law(
    umi_data,
    U_col = "U",
    mi_col = "I_empirical",
    start = umi_start,
    control = umi_control
  )
  joint_fit <- scscale_fit_empirical_joint_law(
    grid,
    n_col = "n",
    U_col = "U",
    mi_col = "I_empirical",
    cell_fit = cell_fit,
    umi_fit = umi_fit,
    mode = joint_mode,
    start = joint_start,
    control = joint_control
  )

  ceiling <- scscale_empirical_global_ceiling(cell_fit, umi_fit, n0 = n_ref, U0 = U_ref)
  joint_coef <- c(joint_fit$coefficients, joint_fit$fixed)
  ceiling <- c(ceiling, joint_fit = joint_coef[["I_infinity"]])
  ceiling <- c(ceiling, max_abs_gap_to_joint = max(abs(ceiling[c("from_cell_marginal", "from_umi_marginal")] - ceiling[["joint_fit"]])))

  fitted <- stats::predict(joint_fit, grid)
  list(
    grid = grid,
    cell_fit = cell_fit,
    umi_fit = umi_fit,
    joint_fit = joint_fit,
    ceiling = ceiling,
    reference = list(n_ref = n_ref, U_ref = U_ref),
    fitted = fitted,
    residuals = grid$I_empirical - fitted,
    metrics = joint_fit$metrics
  )
}

scscale_empirical_scaling_fit <- function(
  x,
  target,
  batch = NULL,
  n_grid = NULL,
  cell_sampling_rates = c(0.125, 0.1875, 0.25, 0.375, 0.50, 0.75, 1.00),
  U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85, 1.00),
  r = 10,
  n_replicates = 10,
  cell_replicates = n_replicates,
  umi_replicates = 1,
  n_workers = 1,
  n_ref = NULL,
  U_ref = NULL,
  joint_mode = c("hybrid", "transfer", "free"),
  correction = c("none", "scvi", "custom"),
  correction_fun = NULL,
  correction_cache = NULL,
  corrected_input = c("normalized", "counts"),
  correction_options = list(),
  seed = 1,
  input = c("counts", "normalized"),
  target_input = input,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  eps = 1e-12,
  use_irlba = TRUE,
  cell_start = NULL,
  umi_start = NULL,
  joint_start = NULL,
  cell_control = list(maxit = 10000),
  umi_control = list(maxit = 10000),
  joint_control = list(maxit = 20000)
) {
  joint_mode <- match.arg(joint_mode)
  correction <- match.arg(correction)
  input <- match.arg(input)
  corrected_input <- match.arg(corrected_input)
  target_input <- match.arg(target_input, choices = c("counts", "normalized"))
  count_transform <- match.arg(count_transform)
  if (!is.null(U_grid) && missing(sampling_rates)) {
    sampling_rates <- NULL
  }
  if (!is.null(n_grid) && missing(cell_sampling_rates)) {
    cell_sampling_rates <- NULL
  }

  if (correction == "none") {
    grid <- scscale_empirical_inu_grid(
      x = x,
      target = target,
      n_grid = n_grid,
      cell_sampling_rates = cell_sampling_rates,
      U_grid = U_grid,
      sampling_rates = sampling_rates,
      r = r,
      n_replicates = cell_replicates,
      n_workers = n_workers,
      seed = seed,
      input = input,
      target_input = target_input,
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale,
      eps = eps,
      use_irlba = use_irlba
    )
    if (!"cell_replicate" %in% names(grid)) grid$cell_replicate <- grid$replicate
    if (!"umi_replicate" %in% names(grid)) grid$umi_replicate <- 1L
    grid$correction <- "none"
  } else {
    grid <- scscale_empirical_inu_grid_corrected(
      x = x,
      target = target,
      batch = batch,
      n_grid = n_grid,
      cell_sampling_rates = cell_sampling_rates,
      U_grid = U_grid,
      sampling_rates = sampling_rates,
      r = r,
      cell_replicates = cell_replicates,
      umi_replicates = umi_replicates,
      n_workers = n_workers,
      correction = correction,
      correction_fun = correction_fun,
      correction_cache = correction_cache,
      correction_options = correction_options,
      seed = seed,
      input = input,
      corrected_input = corrected_input,
      target_input = target_input,
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale,
      eps = eps,
      use_irlba = use_irlba
    )
  }

  fit <- scscale_fit_empirical_laws_from_grid(
    grid,
    n_ref = n_ref,
    U_ref = U_ref,
    joint_mode = joint_mode,
    cell_start = cell_start,
    umi_start = umi_start,
    joint_start = joint_start,
    cell_control = cell_control,
    umi_control = umi_control,
    joint_control = joint_control
  )

  n_total <- max(grid$n / grid$cell_sampling_rate, na.rm = TRUE)
  U_total <- if (any(is.finite(grid$sampling_rate))) {
    max(grid$U / grid$sampling_rate, na.rm = TRUE)
  } else {
    NA_real_
  }
  out <- list(
    call = match.call(),
    grid = grid,
    cell_fit = fit$cell_fit,
    umi_fit = fit$umi_fit,
    joint_fit = fit$joint_fit,
    ceiling = fit$ceiling,
    reference = fit$reference,
    sampling = list(n_total = n_total, U_total = U_total),
    r = r,
    n_replicates = cell_replicates,
    cell_replicates = cell_replicates,
    umi_replicates = umi_replicates,
    n_workers = n_workers,
    joint_mode = joint_mode,
    correction = list(
      method = correction,
      cache = correction_cache,
      corrected_input = corrected_input,
      options = correction_options
    ),
    fitted = fit$fitted,
    residuals = fit$residuals,
    metrics = fit$metrics,
    preprocessing = list(
      input = input,
      target_input = target_input,
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale
    )
  )
  class(out) <- "scscale_empirical_inu_fit"
  out
}

scscale_empirical_inu_grid <- function(
  x,
  target,
  n_grid = NULL,
  cell_sampling_rates = c(0.125, 0.1875, 0.25, 0.375, 0.50, 0.75, 1.00),
  U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85, 1.00),
  r = 10,
  n_replicates = 10,
  n_workers = 1,
  seed = 1,
  input = c("counts", "normalized"),
  target_input = input,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  eps = 1e-12,
  use_irlba = TRUE
) {
  input <- match.arg(input)
  target_input <- match.arg(target_input, choices = c("counts", "normalized"))
  count_transform <- match.arg(count_transform)
  if (!is.null(U_grid) && missing(sampling_rates)) {
    sampling_rates <- NULL
  }
  if (!is.null(n_grid) && missing(cell_sampling_rates)) {
    cell_sampling_rates <- NULL
  }
  if (is.null(n_grid) && is.null(cell_sampling_rates)) stop("Provide n_grid or cell_sampling_rates.", call. = FALSE)
  if (!is.null(n_grid) && !is.null(cell_sampling_rates)) stop("Provide only one of n_grid or cell_sampling_rates.", call. = FALSE)
  if (is.null(U_grid) && is.null(sampling_rates)) stop("Provide U_grid or sampling_rates.", call. = FALSE)
  if (!is.null(U_grid) && !is.null(sampling_rates)) stop("Provide only one of U_grid or sampling_rates.", call. = FALSE)

  x <- counts_matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
  n_replicates <- as.integer(n_replicates)
  if (!is.finite(n_replicates) || n_replicates < 1L) stop("n_replicates must be positive.", call. = FALSE)
  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) stop("n_workers must be positive.", call. = FALSE)

  target_is_matrix <- is.matrix(target) || inherits(target, "Matrix")
  if (target_is_matrix) {
    target <- counts_matrix(target)
    if (is.null(colnames(target))) colnames(target) <- colnames(x)
    common <- intersect(colnames(x), colnames(target))
    if (length(common) < 2L) stop("x and matrix target must share at least two cells.", call. = FALSE)
    x <- x[, common, drop = FALSE]
    target <- target[, common, drop = FALSE]
  } else {
    target <- align_vector_to_cells(target, colnames(x), "target")
  }

  if (!is.null(cell_sampling_rates)) {
    cell_sampling_rates <- sort(unique(as.numeric(cell_sampling_rates)))
    if (any(!is.finite(cell_sampling_rates) | cell_sampling_rates <= 0 | cell_sampling_rates > 1)) {
      stop("cell_sampling_rates must be in (0, 1].", call. = FALSE)
    }
    n_grid <- sort(unique(pmax(2L, as.integer(round(cell_sampling_rates * ncol(x))))))
    cell_sampling_rates <- n_grid / ncol(x)
  } else {
    n_grid <- sort(unique(as.integer(n_grid)))
    n_grid <- n_grid[is.finite(n_grid) & n_grid > 1L]
    if (!length(n_grid)) stop("n_grid must contain positive cell counts.", call. = FALSE)
    if (any(n_grid > ncol(x))) stop("n_grid cannot exceed the number of aligned cells.", call. = FALSE)
    cell_sampling_rates <- n_grid / ncol(x)
  }

  if (!is.null(sampling_rates)) {
    sampling_rates <- as.numeric(sampling_rates)
    if (any(!is.finite(sampling_rates) | sampling_rates <= 0 | sampling_rates > 1)) {
      stop("sampling_rates must be in (0, 1].", call. = FALSE)
    }
    U_grid <- sampling_rates * stats::median(colSums(x))
  } else {
    U_grid <- as.numeric(U_grid)
    if (any(!is.finite(U_grid) | U_grid <= 0)) stop("U_grid must contain positive values.", call. = FALSE)
    sampling_rates <- rep(NA_real_, length(U_grid))
  }

  total_umi <- colSums(x)
  eligible_by_u <- vector("list", length(U_grid))
  tasks <- list()
  task_i <- 0L
  for (u_i in seq_along(U_grid)) {
    U <- U_grid[[u_i]]
    rate <- sampling_rates[[u_i]]
    eligible <- if (input == "counts" && is.na(rate)) colnames(x)[total_umi >= U] else colnames(x)
    eligible_by_u[[u_i]] <- eligible
    for (n_target in n_grid) {
      if (length(eligible) < n_target) next
      for (replicate in seq_len(n_replicates)) {
        task_i <- task_i + 1L
        tasks[[task_i]] <- list(u_i = u_i, U = U, rate = rate, n_target = n_target, replicate = replicate)
      }
    }
  }
  if (!length(tasks)) stop("No grid rows were produced; check n_grid and U_grid eligibility.", call. = FALSE)

  compute_row <- function(task) {
    view_seed <- seed + task$u_i * 1000003L + task$n_target * 9176L + task$replicate
    set.seed(view_seed)
    cells <- sample(eligible_by_u[[task$u_i]], task$n_target)
    x_view <- x[, cells, drop = FALSE]
    if (input == "counts") {
      x_view <- if (is.na(task$rate)) {
        scscale_downsample_counts(x_view, U = task$U, seed = view_seed + 13L)
      } else {
        scscale_downsample_counts_fraction(x_view, fraction = task$rate, seed = view_seed + 13L)
      }
      keep <- colSums(x_view) > 0
      x_view <- x_view[, keep, drop = FALSE]
      cells <- colnames(x_view)
    }
    target_view <- if (target_is_matrix) target[, cells, drop = FALSE] else target[cells]
    mi <- scscale_empirical_mi(
      x_view,
      target_view,
      r = r,
      input = input,
      target_input = target_input,
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale,
      eps = eps,
      use_irlba = use_irlba,
      store_subspaces = FALSE
    )
    data.frame(
      n = ncol(x_view),
      cell_sampling_rate = task$n_target / ncol(x),
      U = task$U,
      sampling_rate = task$rate,
      replicate = task$replicate,
      cell_replicate = task$replicate,
      umi_replicate = 1L,
      I_empirical = mi$mi,
      r = r,
      r_eff = if (is.null(mi$r_eff)) NA_integer_ else mi$r_eff,
      actual_median_umi_per_cell = if (input == "counts") stats::median(colSums(x_view)) else NA_real_,
      actual_mean_umi_per_cell = if (input == "counts") mean(colSums(x_view)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  rows <- if (n_workers == 1L || length(tasks) == 1L) {
    lapply(tasks, compute_row)
  } else if (.Platform$OS.type == "unix") {
    parallel::mclapply(tasks, compute_row, mc.cores = min(n_workers, length(tasks)), mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(min(n_workers, length(tasks)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, tasks, compute_row)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) stop("No grid rows were produced; check n_grid and U_grid eligibility.", call. = FALSE)
  out
}

scscale_empirical_inu_grid_corrected <- function(
  x,
  target,
  batch = NULL,
  n_grid = NULL,
  cell_sampling_rates = c(0.125, 0.1875, 0.25, 0.375, 0.50, 0.75, 1.00),
  U_grid = NULL,
  sampling_rates = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.85, 1.00),
  r = 10,
  cell_replicates = 10,
  umi_replicates = 1,
  n_workers = 1,
  correction = c("custom", "scvi"),
  correction_fun = NULL,
  correction_cache = NULL,
  correction_options = list(),
  seed = 1,
  input = c("counts", "normalized"),
  corrected_input = c("normalized", "counts"),
  target_input = input,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  eps = 1e-12,
  use_irlba = TRUE
) {
  correction <- match.arg(correction)
  input <- match.arg(input)
  corrected_input <- match.arg(corrected_input)
  target_input <- match.arg(target_input, choices = c("counts", "normalized"))
  count_transform <- match.arg(count_transform)
  if (correction == "scvi" && is.null(correction_fun)) {
    stop(
      "correction = 'scvi' needs a correction_fun in this build. ",
      "Use a function that runs scVI externally and returns a feature-by-cell corrected matrix.",
      call. = FALSE
    )
  }
  if (!is.function(correction_fun)) {
    stop("correction_fun must be a function for corrected empirical scaling.", call. = FALSE)
  }
  if (input != "counts") {
    stop("corrected empirical scaling requires raw counts input so UMI thinning is well-defined.", call. = FALSE)
  }
  if (!is.null(U_grid) && missing(sampling_rates)) {
    sampling_rates <- NULL
  }
  if (!is.null(n_grid) && missing(cell_sampling_rates)) {
    cell_sampling_rates <- NULL
  }
  if (is.null(n_grid) && is.null(cell_sampling_rates)) stop("Provide n_grid or cell_sampling_rates.", call. = FALSE)
  if (!is.null(n_grid) && !is.null(cell_sampling_rates)) stop("Provide only one of n_grid or cell_sampling_rates.", call. = FALSE)
  if (is.null(U_grid) && is.null(sampling_rates)) stop("Provide U_grid or sampling_rates.", call. = FALSE)
  if (!is.null(U_grid) && !is.null(sampling_rates)) stop("Provide only one of U_grid or sampling_rates.", call. = FALSE)

  x <- counts_matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
  cell_replicates <- as.integer(cell_replicates)
  umi_replicates <- as.integer(umi_replicates)
  n_workers <- as.integer(n_workers)
  if (!is.finite(cell_replicates) || cell_replicates < 1L) stop("cell_replicates must be positive.", call. = FALSE)
  if (!is.finite(umi_replicates) || umi_replicates < 1L) stop("umi_replicates must be positive.", call. = FALSE)
  if (!is.finite(n_workers) || n_workers < 1L) stop("n_workers must be positive.", call. = FALSE)

  target_is_matrix <- is.matrix(target) || inherits(target, "Matrix")
  if (target_is_matrix) {
    target <- counts_matrix(target)
    if (is.null(colnames(target))) colnames(target) <- colnames(x)
    common <- intersect(colnames(x), colnames(target))
    if (length(common) < 2L) stop("x and matrix target must share at least two cells.", call. = FALSE)
    x <- x[, common, drop = FALSE]
    target <- target[, common, drop = FALSE]
  } else {
    target <- align_vector_to_cells(target, colnames(x), "target")
  }
  if (!is.null(batch)) {
    batch <- align_vector_to_cells(batch, colnames(x), "batch")
  }

  if (!is.null(cell_sampling_rates)) {
    cell_sampling_rates <- sort(unique(as.numeric(cell_sampling_rates)))
    if (any(!is.finite(cell_sampling_rates) | cell_sampling_rates <= 0 | cell_sampling_rates > 1)) {
      stop("cell_sampling_rates must be in (0, 1].", call. = FALSE)
    }
    n_grid <- sort(unique(pmax(2L, as.integer(round(cell_sampling_rates * ncol(x))))))
    cell_sampling_rates <- n_grid / ncol(x)
  } else {
    n_grid <- sort(unique(as.integer(n_grid)))
    n_grid <- n_grid[is.finite(n_grid) & n_grid > 1L]
    if (!length(n_grid)) stop("n_grid must contain positive cell counts.", call. = FALSE)
    if (any(n_grid > ncol(x))) stop("n_grid cannot exceed the number of aligned cells.", call. = FALSE)
    cell_sampling_rates <- n_grid / ncol(x)
  }

  if (!is.null(sampling_rates)) {
    sampling_rates <- as.numeric(sampling_rates)
    if (any(!is.finite(sampling_rates) | sampling_rates <= 0 | sampling_rates > 1)) {
      stop("sampling_rates must be in (0, 1].", call. = FALSE)
    }
    U_grid <- sampling_rates * stats::median(colSums(x))
  } else {
    U_grid <- as.numeric(U_grid)
    if (any(!is.finite(U_grid) | U_grid <= 0)) stop("U_grid must contain positive values.", call. = FALSE)
    sampling_rates <- rep(NA_real_, length(U_grid))
  }

  if (!is.null(correction_cache)) {
    dir.create(correction_cache, recursive = TRUE, showWarnings = FALSE)
  }

  all_rows <- list()
  row_i <- 0L
  total_umi <- colSums(x)
  for (u_i in seq_along(U_grid)) {
    U <- U_grid[[u_i]]
    rate <- sampling_rates[[u_i]]
    for (umi_replicate in seq_len(umi_replicates)) {
      view_seed <- seed + u_i * 1000003L + umi_replicate * 65537L
      cache_file <- if (is.null(correction_cache)) {
        NULL
      } else {
        file.path(correction_cache, paste0("corrected_u", u_i, "_rep", umi_replicate, ".rds"))
      }

      if (!is.null(cache_file) && file.exists(cache_file)) {
        corrected <- readRDS(cache_file)
      } else {
        x_u <- if (is.na(rate)) {
          scscale_downsample_counts(x, U = U, seed = view_seed + 13L)
        } else {
          scscale_downsample_counts_fraction(x, fraction = rate, seed = view_seed + 13L)
        }
        keep_u <- colSums(x_u) > 0
        x_u <- x_u[, keep_u, drop = FALSE]
        target_u <- if (target_is_matrix) target[, colnames(x_u), drop = FALSE] else target[colnames(x_u)]
        batch_u <- if (is.null(batch)) NULL else batch[colnames(x_u)]
        corrected <- correction_fun(
          counts = x_u,
          target = target_u,
          batch = batch_u,
          U = U,
          sampling_rate = rate,
          umi_replicate = umi_replicate,
          seed = view_seed,
          options = correction_options
        )
        corrected <- counts_matrix(corrected)
        if (is.null(colnames(corrected))) colnames(corrected) <- colnames(x_u)
        common <- intersect(colnames(corrected), colnames(x_u))
        if (length(common) < 2L) {
          stop("correction_fun must return a feature-by-cell matrix with at least two corrected cells.", call. = FALSE)
        }
        corrected <- corrected[, common, drop = FALSE]
        x_u <- x_u[, common, drop = FALSE]
        if (!is.null(cache_file)) {
          saveRDS(list(corrected = corrected, counts = x_u), cache_file)
        } else {
          corrected <- list(corrected = corrected, counts = x_u)
        }
      }
      if (is.list(corrected) && !is.null(corrected$corrected)) {
        corrected_counts <- corrected$counts
        corrected <- corrected$corrected
      } else {
        corrected_counts <- NULL
      }

      eligible <- colnames(corrected)
      for (n_target in n_grid) {
        if (length(eligible) < n_target) next
        tasks <- lapply(seq_len(cell_replicates), function(cell_replicate) {
          list(n_target = n_target, cell_replicate = cell_replicate)
        })
        compute_row <- function(task) {
          cell_seed <- view_seed + task$n_target * 9176L + task$cell_replicate
          set.seed(cell_seed)
          cells <- sample(eligible, task$n_target)
          x_view <- corrected[, cells, drop = FALSE]
          target_view <- if (target_is_matrix) target[, cells, drop = FALSE] else target[cells]
          mi <- scscale_empirical_mi(
            x_view,
            target_view,
            r = r,
            input = corrected_input,
            target_input = if (target_is_matrix) target_input else corrected_input,
            target_depth = target_depth,
            count_transform = count_transform,
            center = center,
            scale = scale,
            eps = eps,
            use_irlba = use_irlba,
            store_subspaces = FALSE
          )
          raw_view <- if (is.null(corrected_counts)) NULL else corrected_counts[, cells, drop = FALSE]
          data.frame(
            n = ncol(x_view),
            cell_sampling_rate = task$n_target / ncol(x),
            U = U,
            sampling_rate = rate,
            replicate = paste(umi_replicate, task$cell_replicate, sep = "."),
            cell_replicate = task$cell_replicate,
            umi_replicate = umi_replicate,
            correction = correction,
            I_empirical = mi$mi,
            r = r,
            r_eff = if (is.null(mi$r_eff)) NA_integer_ else mi$r_eff,
            actual_median_umi_per_cell = if (is.null(raw_view)) NA_real_ else stats::median(colSums(raw_view)),
            actual_mean_umi_per_cell = if (is.null(raw_view)) NA_real_ else mean(colSums(raw_view)),
            stringsAsFactors = FALSE
          )
        }
        rows <- if (n_workers == 1L || length(tasks) == 1L) {
          lapply(tasks, compute_row)
        } else if (.Platform$OS.type == "unix") {
          parallel::mclapply(tasks, compute_row, mc.cores = min(n_workers, length(tasks)), mc.preschedule = FALSE)
        } else {
          cl <- parallel::makeCluster(min(n_workers, length(tasks)))
          on.exit(parallel::stopCluster(cl), add = TRUE)
          parallel::parLapply(cl, tasks, compute_row)
        }
        for (row in rows) {
          row_i <- row_i + 1L
          all_rows[[row_i]] <- row
        }
      }
    }
  }

  out <- do.call(rbind, all_rows)
  if (is.null(out)) stop("No grid rows were produced; check n_grid and U_grid eligibility.", call. = FALSE)
  out
}

scscale_check_empirical_inu_fit <- function(object) {
  if (!inherits(object, "scscale_empirical_inu_fit")) {
    stop("object must be a scscale_empirical_inu_fit.", call. = FALSE)
  }
  invisible(object)
}

scscale_empirical_n_total <- function(object) {
  n_total <- object$sampling$n_total
  if (is.null(n_total) || !is.finite(n_total)) {
    n_total <- max(object$grid$n / object$grid$cell_sampling_rate, na.rm = TRUE)
  }
  n_total
}

scscale_empirical_U_total <- function(object) {
  U_total <- object$sampling$U_total
  if (is.null(U_total) || !is.finite(U_total)) {
    ok <- is.finite(object$grid$sampling_rate) & object$grid$sampling_rate > 0
    U_total <- if (any(ok)) max(object$grid$U[ok] / object$grid$sampling_rate[ok], na.rm = TRUE) else max(object$grid$U, na.rm = TRUE)
  }
  U_total
}

scscale_resolve_cell_input <- function(object, cell_sampling_rate = NULL, n = NULL, n_points = NULL) {
  n_total <- scscale_empirical_n_total(object)
  if (!is.null(cell_sampling_rate) && !is.null(n)) {
    stop("Provide only one of cell_sampling_rate or n.", call. = FALSE)
  }
  if (is.null(cell_sampling_rate) && is.null(n)) {
    if (is.null(n_points)) {
      cell_sampling_rate <- sort(unique(object$grid$cell_sampling_rate))
    } else {
      cell_sampling_rate <- seq(min(object$grid$cell_sampling_rate), max(object$grid$cell_sampling_rate), length.out = n_points)
    }
    n <- cell_sampling_rate * n_total
  } else if (!is.null(cell_sampling_rate)) {
    cell_sampling_rate <- as.numeric(cell_sampling_rate)
    n <- cell_sampling_rate * n_total
  } else {
    n <- as.numeric(n)
    cell_sampling_rate <- n / n_total
  }
  if (any(!is.finite(n) | n <= 0) || any(!is.finite(cell_sampling_rate) | cell_sampling_rate <= 0)) {
    stop("cell_sampling_rate and n must be positive and finite.", call. = FALSE)
  }
  data.frame(cell_sampling_rate = cell_sampling_rate, n = n)
}

scscale_resolve_umi_input <- function(object, sampling_rate = NULL, U = NULL, n_points = NULL) {
  U_total <- scscale_empirical_U_total(object)
  if (!is.null(sampling_rate) && !is.null(U)) {
    stop("Provide only one of sampling_rate or U.", call. = FALSE)
  }
  if (is.null(sampling_rate) && is.null(U)) {
    rates <- object$grid$sampling_rate
    rates <- rates[is.finite(rates)]
    if (is.null(n_points)) {
      sampling_rate <- sort(unique(rates))
    } else {
      sampling_rate <- seq(min(rates), max(rates), length.out = n_points)
    }
    U <- sampling_rate * U_total
  } else if (!is.null(sampling_rate)) {
    sampling_rate <- as.numeric(sampling_rate)
    U <- sampling_rate * U_total
  } else {
    U <- as.numeric(U)
    sampling_rate <- U / U_total
  }
  if (any(!is.finite(U) | U <= 0) || any(!is.finite(sampling_rate) | sampling_rate <= 0)) {
    stop("sampling_rate and U must be positive and finite.", call. = FALSE)
  }
  data.frame(sampling_rate = sampling_rate, U = U)
}

scscale_empirical_summary <- function(data, by, y = "I_empirical") {
  mean_df <- stats::aggregate(stats::as.formula(paste(y, "~", paste(by, collapse = " + "))), data = data, FUN = mean)
  sd_df <- stats::aggregate(stats::as.formula(paste(y, "~", paste(by, collapse = " + "))), data = data, FUN = stats::sd)
  n_df <- stats::aggregate(stats::as.formula(paste("replicate ~", paste(by, collapse = " + "))), data = data, FUN = function(x) length(unique(x)))
  mean_df$I_sd <- sd_df[[y]]
  mean_df$n_replicates <- n_df$replicate
  mean_df
}

scscale_draw_interval_bars <- function(x, y, sd, col, cap = NULL, lwd = 1.6) {
  sd[!is.finite(sd)] <- 0
  if (is.null(cap)) cap <- diff(range(x)) * 0.012
  graphics::segments(x, y - sd, x, y + sd, col = col, lwd = lwd)
  graphics::segments(x - cap, y - sd, x + cap, y - sd, col = col, lwd = lwd)
  graphics::segments(x - cap, y + sd, x + cap, y + sd, col = col, lwd = lwd)
}

scscale_empirical_values <- function(object, summarize = FALSE) {
  scscale_check_empirical_inu_fit(object)
  out <- object$grid
  out$I_fit <- object$fitted
  out$residual <- object$residuals
  if (!summarize) return(out)
  summary <- scscale_empirical_summary(out, by = c("cell_sampling_rate", "sampling_rate", "U"))
  fit_summary <- stats::aggregate(I_fit ~ cell_sampling_rate + sampling_rate + U, data = out, FUN = mean)
  residual_summary <- stats::aggregate(residual ~ cell_sampling_rate + sampling_rate + U, data = out, FUN = mean)
  summary$I_fit <- fit_summary$I_fit
  summary$residual <- residual_summary$residual
  summary
}

scscale_empirical_parameters <- function(object) {
  scscale_check_empirical_inu_fit(object)
  rows <- list(
    cell = object$cell_fit$coefficients,
    umi = object$umi_fit$coefficients,
    joint = c(object$joint_fit$coefficients, object$joint_fit$fixed),
    ceiling = object$ceiling
  )
  out <- do.call(rbind, lapply(names(rows), function(model) {
    data.frame(model = model, parameter = names(rows[[model]]), value = as.numeric(rows[[model]]), row.names = NULL)
  }))
  row.names(out) <- NULL
  out
}

scscale_predict_empirical_cell <- function(object, cell_sampling_rate = NULL, n = NULL, n_points = NULL) {
  scscale_check_empirical_inu_fit(object)
  out <- scscale_resolve_cell_input(object, cell_sampling_rate = cell_sampling_rate, n = n, n_points = n_points)
  out$I_fit <- stats::predict(object$cell_fit, out)
  out
}

scscale_predict_empirical_umi <- function(object, sampling_rate = NULL, U = NULL, n_points = NULL) {
  scscale_check_empirical_inu_fit(object)
  out <- scscale_resolve_umi_input(object, sampling_rate = sampling_rate, U = U, n_points = n_points)
  out$I_fit <- stats::predict(object$umi_fit, out)
  out
}

scscale_predict_empirical_joint <- function(
  object,
  cell_sampling_rate = NULL,
  n = NULL,
  sampling_rate = NULL,
  U = NULL,
  cell_points = NULL,
  umi_points = NULL
) {
  scscale_check_empirical_inu_fit(object)
  cell <- scscale_resolve_cell_input(object, cell_sampling_rate = cell_sampling_rate, n = n, n_points = cell_points)
  umi <- scscale_resolve_umi_input(object, sampling_rate = sampling_rate, U = U, n_points = umi_points)
  out <- merge(cell, umi, by = NULL)
  out$I_fit <- stats::predict(object, out)
  out
}

scscale_plot_empirical_fit <- function(object, type = c("both", "cell", "umi"), n_points = 200, ...) {
  scscale_check_empirical_inu_fit(object)
  type <- match.arg(type)
  old_par <- NULL
  if (type == "both") {
    old_par <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.3, 2.8, 1))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  if (type %in% c("both", "cell")) {
    obs <- object$grid[abs(object$grid$U - object$reference$U_ref) < 1e-8, , drop = FALSE]
    obs <- scscale_empirical_summary(obs, by = "cell_sampling_rate")
    curve <- scscale_predict_empirical_cell(object, n_points = n_points)
    I_infinity <- c(object$cell_fit$coefficients, object$cell_fit$fixed)[["I_infinity"]]
    graphics::plot(
      NA,
      xlim = range(obs$cell_sampling_rate),
      ylim = range(c(obs$I_empirical - obs$I_sd, obs$I_empirical + obs$I_sd, curve$I_fit, I_infinity), na.rm = TRUE),
      xlab = "Cell sampling rate",
      ylab = "Empirical MI",
      main = "Cell-number scaling",
      ...
    )
    scscale_draw_interval_bars(obs$cell_sampling_rate, obs$I_empirical, obs$I_sd, col = "#0072B2")
    graphics::lines(curve$cell_sampling_rate, curve$I_fit, col = "#D55E00", lwd = 2)
    graphics::abline(h = I_infinity, col = "grey35", lwd = 2, lty = 2)
    graphics::legend(
      "bottomright",
      legend = expression("empirical mean +/- sd", "curve fit", I[infinity] ~ "upper bound"),
      col = c("#0072B2", "#D55E00", "grey35"),
      lwd = c(1.6, 2, 2),
      lty = c(1, 1, 2),
      bty = "n"
    )
  }
  if (type %in% c("both", "umi")) {
    obs <- object$grid[object$grid$n == object$reference$n_ref, , drop = FALSE]
    obs <- scscale_empirical_summary(obs, by = c("sampling_rate", "U"))
    curve <- scscale_predict_empirical_umi(object, n_points = n_points)
    I_infinity <- c(object$umi_fit$coefficients, object$umi_fit$fixed)[["I_infinity"]]
    graphics::plot(
      NA,
      xlim = range(obs$sampling_rate),
      ylim = range(c(obs$I_empirical - obs$I_sd, obs$I_empirical + obs$I_sd, curve$I_fit, I_infinity), na.rm = TRUE),
      xlab = "UMI sampling rate",
      ylab = "Empirical MI",
      main = "UMI scaling",
      ...
    )
    scscale_draw_interval_bars(obs$sampling_rate, obs$I_empirical, obs$I_sd, col = "#0072B2")
    graphics::lines(curve$sampling_rate, curve$I_fit, col = "#D55E00", lwd = 2)
    graphics::abline(h = I_infinity, col = "grey35", lwd = 2, lty = 2)
    graphics::legend(
      "bottomright",
      legend = expression("empirical mean +/- sd", "curve fit", I[infinity] ~ "upper bound"),
      col = c("#0072B2", "#D55E00", "grey35"),
      lwd = c(1.6, 2, 2),
      lty = c(1, 1, 2),
      bty = "n"
    )
  }
  invisible(object)
}

scscale_plot_empirical_umi_levels <- function(object, sampling_rate = NULL, U = NULL, n_points = 300, ...) {
  scscale_check_empirical_inu_fit(object)
  if (is.null(sampling_rate) && is.null(U)) sampling_rate <- c(0.20, 0.50, 1.00)
  umi <- scscale_resolve_umi_input(object, sampling_rate = sampling_rate, U = U)
  obs <- scscale_empirical_summary(object$grid, by = c("cell_sampling_rate", "sampling_rate", "U"))
  cols <- grDevices::hcl.colors(nrow(umi), palette = "Dark 3")
  curve <- scscale_predict_empirical_joint(object, sampling_rate = umi$sampling_rate, cell_points = n_points)
  graphics::plot(
    NA,
    xlim = range(object$grid$cell_sampling_rate),
    ylim = range(c(obs$I_empirical - obs$I_sd, obs$I_empirical + obs$I_sd, curve$I_fit), na.rm = TRUE),
    xlab = "Cell sampling rate",
    ylab = "Empirical MI",
    main = "Joint scaling by UMI level",
    ...
  )
  for (i in seq_len(nrow(umi))) {
    obs_i <- obs[abs(obs$U - umi$U[[i]]) < 1e-8, , drop = FALSE]
    obs_i <- obs_i[order(obs_i$cell_sampling_rate), ]
    if (nrow(obs_i)) {
      scscale_draw_interval_bars(obs_i$cell_sampling_rate, obs_i$I_empirical, obs_i$I_sd, col = cols[[i]])
    }
    curve_i <- curve[abs(curve$U - umi$U[[i]]) < 1e-8, , drop = FALSE]
    curve_i <- curve_i[order(curve_i$cell_sampling_rate), ]
    graphics::lines(curve_i$cell_sampling_rate, curve_i$I_fit, col = cols[[i]], lwd = 2)
  }
  graphics::legend("bottomright", legend = paste0("rho = ", signif(umi$sampling_rate, 3)), col = cols, lwd = 2, bty = "n", cex = 0.85)
  invisible(object)
}

scscale_plot_empirical_surface <- function(
  object,
  n = NULL,
  U = NULL,
  n_points = 60,
  U_points = 60,
  show_points = TRUE,
  theta = 35,
  phi = 25,
  expand = 0.65,
  surface_col = "#56B4E9",
  point_col = "#D55E00",
  border = NA,
  shade = 0.35,
  ...
) {
  scscale_check_empirical_inu_fit(object)
  if (is.null(n)) {
    n <- seq(min(object$grid$n), max(object$grid$n), length.out = n_points)
  }
  if (is.null(U)) {
    U <- seq(min(object$grid$U), max(object$grid$U), length.out = U_points)
  }
  n <- sort(unique(as.numeric(n)))
  U <- sort(unique(as.numeric(U)))
  if (length(n) < 2L || any(!is.finite(n) | n <= 0)) {
    stop("n must contain at least two unique positive finite cell numbers.", call. = FALSE)
  }
  if (length(U) < 2L || any(!is.finite(U) | U <= 0)) {
    stop("U must contain at least two unique positive finite UMI depths.", call. = FALSE)
  }

  surface <- expand.grid(n = n, U = U)
  surface$I_fit <- stats::predict(object, surface)
  z <- matrix(surface$I_fit, nrow = length(n), ncol = length(U))
  z_values <- z
  observed <- NULL
  if (isTRUE(show_points)) {
    observed <- scscale_empirical_summary(object$grid, by = c("n", "U"))
    observed <- observed[
      observed$n >= min(n) & observed$n <= max(n) &
        observed$U >= min(U) & observed$U <= max(U),
      ,
      drop = FALSE
    ]
    z_values <- c(z_values, observed$I_empirical)
  }

  perspective <- graphics::persp(
    x = n,
    y = U,
    z = z,
    zlim = range(z_values, finite = TRUE),
    theta = theta,
    phi = phi,
    expand = expand,
    col = surface_col,
    border = border,
    shade = shade,
    ticktype = "detailed",
    xlab = "Cell number n",
    ylab = "UMI depth U",
    zlab = "Mutual information I(n,U)",
    main = "Joint scaling surface",
    ...
  )
  if (isTRUE(show_points) && nrow(observed)) {
    projected <- grDevices::trans3d(
      observed$n,
      observed$U,
      observed$I_empirical,
      pmat = perspective
    )
    graphics::points(projected, pch = 19, col = point_col)
  }
  invisible(list(surface = surface, perspective = perspective, observed = observed))
}

predict.scscale_empirical_scaling_fit <- function(object, newdata = NULL, ...) {
  coef <- c(object$coefficients, object$fixed)
  if (is.null(newdata)) newdata <- object$data
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (object$model == "cell") {
    n <- as.numeric(newdata[[object$columns$n]])
    return(scscale_empirical_cell_curve(n, coef[["I_infinity"]], coef[["C_n"]], coef[["A_n"]]))
  }
  if (object$model == "umi") {
    U <- as.numeric(newdata[[object$columns$U]])
    return(scscale_empirical_umi_curve(U, coef[["I_infinity"]], coef[["C_U"]], coef[["A_U"]], coef[["U_offset"]]))
  }
  n <- as.numeric(newdata[[object$columns$n]])
  U <- as.numeric(newdata[[object$columns$U]])
  scscale_empirical_joint_curve(n, U, coef[["I_infinity"]], coef[["C_n"]], coef[["A_n"]], coef[["C_U"]], coef[["A_U"]], coef[["U_offset"]])
}

print.scscale_empirical_scaling_fit <- function(x, ...) {
  cat("scScale empirical scaling fit\n")
  cat("  model: ", x$model, if (!is.null(x$mode)) paste0(" (", x$mode, ")") else "", "\n", sep = "")
  cat("  formula: ", x$formula, "\n", sep = "")
  cat("  coefficients:\n")
  print(signif(x$coefficients, 5))
  if (length(x$fixed)) {
    cat("  fixed:\n")
    print(signif(x$fixed, 5))
  }
  cat("  rmse: ", signif(x$metrics$rmse, 5), "  convergence: ", x$convergence, "\n", sep = "")
  invisible(x)
}

predict.scscale_empirical_inu_fit <- function(object, newdata = NULL, ...) {
  stats::predict(object$joint_fit, newdata = newdata, ...)
}

print.scscale_empirical_inu_fit <- function(x, ...) {
  cat("scScale direct empirical I(n,U) fit\n")
  cat("  grid rows: ", nrow(x$grid), "\n", sep = "")
  cat("  rank: ", x$r, "  cell replicates: ", x$cell_replicates %||% x$n_replicates, "  umi replicates: ", x$umi_replicates %||% 1L, "  workers: ", x$n_workers, "\n", sep = "")
  if (!is.null(x$correction$method)) {
    cat("  correction: ", x$correction$method, "\n", sep = "")
  }
  cat("  reference: n_ref=", x$reference$n_ref, ", U_ref=", signif(x$reference$U_ref, 5), "\n", sep = "")
  cat("  joint mode: ", x$joint_mode, "\n", sep = "")
  coef <- c(x$joint_fit$coefficients, x$joint_fit$fixed)
  cat("  I_infinity: ", signif(coef[["I_infinity"]], 5), "\n", sep = "")
  cat("  joint rmse: ", signif(x$metrics$rmse, 5), "\n", sep = "")
  cat("  ceiling diagnostics:\n")
  print(signif(x$ceiling, 5))
  invisible(x)
}
