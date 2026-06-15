#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

out_dir <- arg_value("--out-dir", "outputs/comparisons/scscale_batch_effect_vs_github")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in sort(list.files("packages/scScale/R", pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

# Oracle helpers copied from the GitHub batch-effect fitting script, with only
# names changed and global settings made explicit for a deterministic smoke test.
gh_build_matrix <- function(X_counts, transform = c("log1p", "clr", "none"), features = NULL) {
  transform <- match.arg(transform)
  if (!is.null(features)) X_counts <- X_counts[intersect(features, rownames(X_counts)), , drop = FALSE]
  X <- as.matrix(X_counts)
  if (identical(transform, "clr")) {
    X_use <- apply(X, 2, function(x) {
      log1p(x / exp(mean(log1p(x))))
    })
  } else if (identical(transform, "log1p")) {
    X_use <- log1p(X)
  } else {
    X_use <- X
  }
  X_cells <- t(X_use)
  X_cells <- scale(X_cells, center = TRUE, scale = TRUE)
  X_cells[!is.finite(X_cells)] <- 0
  out <- t(X_cells)
  rownames(out) <- rownames(X_counts)
  colnames(out) <- colnames(X_counts)
  out
}

gh_get_right_singular_vectors <- function(X, r = 20) {
  X <- as.matrix(X)
  r_eff <- min(r, nrow(X) - 1, ncol(X) - 1)
  if (r_eff < 1) return(NULL)
  if (requireNamespace("irlba", quietly = TRUE) && r_eff < min(dim(X)) - 1L) {
    sv <- irlba::irlba(X, nv = r_eff, nu = 0, maxit = 2000, work = max(3 * r_eff, 30))
    Z <- sv$v[, seq_len(r_eff), drop = FALSE]
  } else {
    sv <- svd(X, nu = 0, nv = r_eff)
    Z <- sv$v[, seq_len(r_eff), drop = FALSE]
  }
  rownames(Z) <- colnames(X)
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  Z
}

gh_orthonormalize_cell_matrix <- function(Z) {
  Z <- as.matrix(Z)
  Z <- scale(Z, center = TRUE, scale = FALSE)
  Z[!is.finite(Z)] <- 0
  keep_col <- apply(Z, 2, function(x) stats::sd(x) > 1e-12)
  Z <- Z[, keep_col, drop = FALSE]
  if (ncol(Z) < 1) return(NULL)
  qrZ <- qr(Z)
  Q <- qr.Q(qrZ)
  Q <- Q[, seq_len(qrZ$rank), drop = FALSE]
  rownames(Q) <- rownames(Z)
  Q
}

gh_compute_subspace_overlap_mi <- function(Z1, Z2, eps = 1e-12) {
  common_cells <- intersect(rownames(Z1), rownames(Z2))
  if (length(common_cells) < 5) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }
  Z1 <- Z1[common_cells, , drop = FALSE]
  Z2 <- Z2[common_cells, , drop = FALSE]
  Q1 <- gh_orthonormalize_cell_matrix(Z1)
  Q2 <- gh_orthonormalize_cell_matrix(Z2)
  if (is.null(Q1) || is.null(Q2)) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }
  r_eff <- min(ncol(Q1), ncol(Q2))
  if (r_eff < 1) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }
  Cmat <- crossprod(Q1, Q2)
  gamma <- svd(Cmat, nu = 0, nv = 0)$d^2
  gamma <- gamma[seq_len(min(r_eff, length(gamma)))]
  gamma <- pmin(pmax(gamma, 0), 1 - eps)
  mi <- -0.5 * sum(log(1 - gamma))
  list(mi = mi, mi_norm = mi / length(gamma), gamma = gamma, r_eff = length(gamma))
}

gh_fit_stats <- function(y, pred) {
  resid <- y - pred
  denom <- sum((y - mean(y))^2)
  data.frame(
    R2 = ifelse(denom > 0, 1 - sum(resid^2) / denom, NA_real_),
    RMSE = sqrt(mean(resid^2)),
    MAE = mean(abs(resid))
  )
}

gh_try_nls <- function(expr, dat, starts, lower, upper) {
  best_fit <- NULL
  best_rmse <- Inf
  errors <- character(0)
  for (ii in seq_along(starts)) {
    fit <- tryCatch(
      stats::nls(
        expr,
        data = dat,
        start = starts[[ii]],
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
        best_fit <- fit
        best_rmse <- rmse
      }
    }
  }
  if (is.null(best_fit)) return(list(fit = NULL, error = paste(unique(errors), collapse = " | ")))
  list(fit = best_fit, error = NA_character_)
}

gh_fit_batch_number_law <- function(summary_df, min_fit_points = 5) {
  dat <- summary_df[is.finite(summary_df$I_fit) & summary_df$m_batch > 0, , drop = FALSE]
  dat <- dat[order(dat$m_batch), , drop = FALSE]
  if (nrow(dat) < min_fit_points) return(NULL)
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
  fit_out <- gh_try_nls(
    I_fit ~ I_inf - C * log(1 - A / m_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = -max_abs_a),
    upper = c(I_inf = Inf, C = Inf, A = max_abs_a)
  )
  if (is.null(fit_out$fit)) return(NULL)
  co <- stats::coef(fit_out$fit)
  dat$I_pred <- as.numeric(stats::predict(fit_out$fit, newdata = dat))
  list(coef = co, data = dat, stats = gh_fit_stats(dat$I_fit, dat$I_pred))
}

gh_fit_cells_per_batch_law <- function(summary_df, min_fit_points = 5) {
  dat <- summary_df[is.finite(summary_df$I_fit) & summary_df$cells_per_batch > 0, , drop = FALSE]
  dat <- dat[order(dat$cells_per_batch), , drop = FALSE]
  if (nrow(dat) < min_fit_points) return(NULL)
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
  fit_out <- gh_try_nls(
    I_fit ~ I_inf - C * log(1 + A / cells_per_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = 0),
    upper = c(I_inf = Inf, C = Inf, A = 20 * s_max)
  )
  if (is.null(fit_out$fit)) return(NULL)
  co <- stats::coef(fit_out$fit)
  dat$I_pred <- as.numeric(stats::predict(fit_out$fit, newdata = dat))
  list(coef = co, data = dat, stats = gh_fit_stats(dat$I_fit, dat$I_pred))
}

max_abs_diff <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(0)
  max(abs(a[ok] - b[ok]))
}

set.seed(20260614)
n_cells <- 48
rna <- matrix(rpois(40 * n_cells, lambda = 4), nrow = 40)
adt <- matrix(rpois(12 * n_cells, lambda = 3), nrow = 12)
latent <- rnorm(n_cells)
rna[1:8, ] <- rna[1:8, ] + rep(latent > 0, each = 8)
adt[1:4, ] <- adt[1:4, ] + rep(latent > 0, each = 4)
colnames(rna) <- colnames(adt) <- paste0("cell_", seq_len(n_cells))
rownames(rna) <- paste0("gene_", seq_len(nrow(rna)))
rownames(adt) <- paste0("adt_", seq_len(nrow(adt)))

features <- rownames(rna)[seq_len(30)]
gh_X <- gh_build_matrix(rna, transform = "log1p", features = features)
gh_Y <- gh_build_matrix(adt, transform = "log1p")
set.seed(99)
gh_Zx <- gh_get_right_singular_vectors(gh_X, r = 5)
gh_Zy <- gh_get_right_singular_vectors(gh_Y, r = 5)
gh_mi <- gh_compute_subspace_overlap_mi(gh_Zx, gh_Zy)

set.seed(99)
pkg_mi <- batch_effect_mi(
  rna[features, , drop = FALSE],
  adt,
  r_x = 5,
  r_y = 5,
  transform_x = "log1p",
  transform_y = "log1p"
)

df_m <- data.frame(
  m_batch = c(2, 3, 5, 8, 12),
  mean_I_bio_norm = c(0.12, 0.15, 0.18, 0.20, 0.215)
)
gh_m_fit <- gh_fit_batch_number_law(transform(df_m, I_fit = mean_I_bio_norm), min_fit_points = 4)
pkg_m_fit <- fit_batch_effect_scaling(df_m, law = "batch_number", target_col = "mean_I_bio_norm", min_points = 4)

df_s <- data.frame(
  cells_per_batch = c(20, 30, 50, 80, 120),
  mean_I_bio_norm = c(0.075, 0.12, 0.17, 0.205, 0.23)
)
gh_s_fit <- gh_fit_cells_per_batch_law(transform(df_s, I_fit = mean_I_bio_norm), min_fit_points = 4)
pkg_s_fit <- fit_batch_effect_scaling(df_s, law = "cells_per_batch", target_col = "mean_I_bio_norm", min_points = 4)

summary <- data.frame(
  diff_mi = abs(gh_mi$mi - pkg_mi$mi),
  diff_mi_norm = abs(gh_mi$mi_norm - pkg_mi$mi_norm),
  max_abs_diff_gamma = max_abs_diff(gh_mi$gamma, pkg_mi$gamma),
  diff_r_eff = gh_mi$r_eff - pkg_mi$r_eff,
  max_abs_diff_batch_number_coef = max_abs_diff(gh_m_fit$coef, pkg_m_fit$coef),
  max_abs_diff_batch_number_pred = max_abs_diff(gh_m_fit$data$I_pred, pkg_m_fit$data$I_pred),
  max_abs_diff_cells_per_batch_coef = max_abs_diff(gh_s_fit$coef, pkg_s_fit$coef),
  max_abs_diff_cells_per_batch_pred = max_abs_diff(gh_s_fit$data$I_pred, pkg_s_fit$data$I_pred)
)

utils::write.csv(summary, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
print(summary)

if (any(unlist(summary) != 0)) {
  stop("scScale batch-effect helpers differ from copied GitHub formulas; see validation outputs.", call. = FALSE)
}

cat("Exact match to copied GitHub batch-effect formulas.\n")
