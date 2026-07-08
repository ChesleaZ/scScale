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

scscale_empirical_scaling_fit <- function(
  x,
  target,
  n_grid,
  U_grid = NULL,
  sampling_rates = NULL,
  r = 10,
  n_replicates = 1,
  n_ref = NULL,
  U_ref = NULL,
  joint_mode = c("hybrid", "transfer", "free"),
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
  input <- match.arg(input)
  target_input <- match.arg(target_input, choices = c("counts", "normalized"))
  count_transform <- match.arg(count_transform)

  grid <- scscale_empirical_inu_grid(
    x = x,
    target = target,
    n_grid = n_grid,
    U_grid = U_grid,
    sampling_rates = sampling_rates,
    r = r,
    n_replicates = n_replicates,
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
  out <- list(
    call = match.call(),
    grid = grid,
    cell_fit = cell_fit,
    umi_fit = umi_fit,
    joint_fit = joint_fit,
    ceiling = ceiling,
    reference = list(n_ref = n_ref, U_ref = U_ref),
    r = r,
    n_replicates = n_replicates,
    joint_mode = joint_mode,
    fitted = fitted,
    residuals = grid$I_empirical - fitted,
    metrics = joint_fit$metrics,
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
  n_grid,
  U_grid = NULL,
  sampling_rates = NULL,
  r = 10,
  n_replicates = 1,
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
  if (is.null(U_grid) && is.null(sampling_rates)) stop("Provide U_grid or sampling_rates.", call. = FALSE)
  if (!is.null(U_grid) && !is.null(sampling_rates)) stop("Provide only one of U_grid or sampling_rates.", call. = FALSE)

  x <- counts_matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
  n_grid <- sort(unique(as.integer(n_grid)))
  n_grid <- n_grid[is.finite(n_grid) & n_grid > 1L]
  if (!length(n_grid)) stop("n_grid must contain positive cell counts.", call. = FALSE)
  n_replicates <- as.integer(n_replicates)
  if (!is.finite(n_replicates) || n_replicates < 1L) stop("n_replicates must be positive.", call. = FALSE)

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
  rows <- list()
  row_i <- 0L
  for (u_i in seq_along(U_grid)) {
    U <- U_grid[[u_i]]
    rate <- sampling_rates[[u_i]]
    eligible <- if (input == "counts" && is.na(rate)) colnames(x)[total_umi >= U] else colnames(x)
    for (n_target in n_grid) {
      if (length(eligible) < n_target) next
      for (replicate in seq_len(n_replicates)) {
        view_seed <- seed + u_i * 1000003L + n_target * 9176L + replicate
        set.seed(view_seed)
        cells <- sample(eligible, n_target)
        x_view <- x[, cells, drop = FALSE]
        if (input == "counts") {
          x_view <- if (is.na(rate)) {
            scscale_downsample_counts(x_view, U = U, seed = view_seed + 13L)
          } else {
            scscale_downsample_counts_fraction(x_view, fraction = rate, seed = view_seed + 13L)
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
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          n = ncol(x_view),
          U = U,
          sampling_rate = rate,
          replicate = replicate,
          I_empirical = mi$mi,
          r = r,
          r_eff = mi$r_eff %||% NA_integer_,
          actual_median_umi_per_cell = if (input == "counts") stats::median(colSums(x_view)) else NA_real_,
          actual_mean_umi_per_cell = if (input == "counts") mean(colSums(x_view)) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) stop("No grid rows were produced; check n_grid and U_grid eligibility.", call. = FALSE)
  out
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
  cat("  rank: ", x$r, "  replicates: ", x$n_replicates, "\n", sep = "")
  cat("  reference: n_ref=", x$reference$n_ref, ", U_ref=", signif(x$reference$U_ref, 5), "\n", sep = "")
  cat("  joint mode: ", x$joint_mode, "\n", sep = "")
  coef <- c(x$joint_fit$coefficients, x$joint_fit$fixed)
  cat("  I_infinity: ", signif(coef[["I_infinity"]], 5), "\n", sep = "")
  cat("  joint rmse: ", signif(x$metrics$rmse, 5), "\n", sep = "")
  cat("  ceiling diagnostics:\n")
  print(signif(x$ceiling, 5))
  invisible(x)
}
