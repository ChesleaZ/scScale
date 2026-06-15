#!/usr/bin/env Rscript

# Spectral overlap mutual information utilities. These functions operate on
# low-rank subspace overlaps, not on discrete cluster labels.

subspace_gamma <- function(left_basis, right_basis, r = 10) {
  r_eff <- min(r, ncol(left_basis), ncol(right_basis))
  if (r_eff < 1) return(numeric())
  s <- t(left_basis[, seq_len(r_eff), drop = FALSE]) %*%
    right_basis[, seq_len(r_eff), drop = FALSE]
  gamma <- svd(s, nu = 0, nv = 0)$d^2
  pmin(pmax(gamma, 0), 1 - 1e-12)
}

overlap_mi_from_gamma <- function(gamma) {
  gamma <- gamma[is.finite(gamma) & gamma > 0]
  if (!length(gamma)) return(0)
  gamma <- pmin(pmax(gamma, 0), 1 - 1e-12)
  -0.5 * sum(log1p(-gamma))
}

cell_subspace_overlap_mi <- function(x_current, x_reference, r = 10) {
  common_genes <- intersect(rownames(x_current), rownames(x_reference))
  common_cells <- intersect(colnames(x_current), colnames(x_reference))
  if (length(common_genes) < 2 || length(common_cells) < 2) {
    return(list(mi = NA_real_, gamma = rep(NA_real_, r), r_eff = 0))
  }

  x_cur <- x_current[common_genes, common_cells, drop = FALSE]
  x_ref <- x_reference[common_genes, common_cells, drop = FALSE]
  r_eff <- min(r, ncol(x_cur), nrow(x_cur), ncol(x_ref), nrow(x_ref))
  if (r_eff < 1) return(list(mi = NA_real_, gamma = rep(NA_real_, r), r_eff = 0))

  v_cur <- svd(x_cur, nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  v_ref <- svd(x_ref, nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  gamma <- subspace_gamma(v_cur, v_ref, r = r_eff)

  gamma_out <- rep(NA_real_, r)
  gamma_out[seq_along(gamma)] <- gamma
  list(mi = overlap_mi_from_gamma(gamma), gamma = gamma_out, r_eff = r_eff)
}

gene_subspace_overlap_mi <- function(x_current, x_reference, r = 10) {
  common_genes <- intersect(rownames(x_current), rownames(x_reference))
  if (length(common_genes) < 2) {
    return(list(mi = NA_real_, gamma = rep(NA_real_, r), r_eff = 0))
  }

  x_cur <- x_current[common_genes, , drop = FALSE]
  x_ref <- x_reference[common_genes, , drop = FALSE]
  r_eff <- min(r, nrow(x_cur) - 1L, ncol(x_cur), nrow(x_ref) - 1L, ncol(x_ref))
  if (r_eff < 1) return(list(mi = NA_real_, gamma = rep(NA_real_, r), r_eff = 0))

  v_cur <- svd(t(x_cur), nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  v_ref <- svd(t(x_ref), nu = 0, nv = r_eff)$v[, seq_len(r_eff), drop = FALSE]
  gamma <- subspace_gamma(v_cur, v_ref, r = r_eff)

  gamma_out <- rep(NA_real_, r)
  gamma_out[seq_along(gamma)] <- gamma
  list(mi = overlap_mi_from_gamma(gamma), gamma = gamma_out, r_eff = r_eff)
}

theoretical_gamma_cells_single <- function(d, gamma, r = 10) {
  d_use <- d[is.finite(d) & d > 0]
  d_use <- d_use[seq_len(min(r, length(d_use)))]
  out <- rep(NA_real_, r)
  if (!length(d_use)) return(out)
  g <- (d_use^4 * gamma - 1) / (d_use^2 * (gamma * d_use^2 + 1))
  g[d_use^4 <= 1 / gamma] <- 0
  out[seq_along(g)] <- pmin(pmax(g, 0), 1 - 1e-12)
  out
}

theoretical_gamma_genes_single <- function(d, gamma, r = 10) {
  d_use <- d[is.finite(d) & d > 0]
  d_use <- d_use[seq_len(min(r, length(d_use)))]
  out <- rep(NA_real_, r)
  if (!length(d_use)) return(out)
  g <- (d_use^4 * gamma - 1) / (gamma * d_use^2 * (d_use^2 + 1))
  g[d_use^4 <= 1 / gamma] <- 0
  out[seq_along(g)] <- pmin(pmax(g, 0), 1 - 1e-12)
  out
}

theoretical_overlap_mi_cells <- function(d_current, d_reference, gamma_current, gamma_reference = gamma_current, r = 10) {
  gamma_cur <- theoretical_gamma_cells_single(d_current, gamma_current, r = r)
  gamma_ref <- theoretical_gamma_cells_single(d_reference, gamma_reference, r = r)
  gamma_cur[!is.finite(gamma_cur)] <- 0
  gamma_ref[!is.finite(gamma_ref)] <- 0
  gamma_double <- pmin(pmax(gamma_cur * gamma_ref, 0), 1 - 1e-12)
  list(
    mi = overlap_mi_from_gamma(gamma_double),
    gamma_current = gamma_cur,
    gamma_reference = gamma_ref,
    gamma_double = gamma_double
  )
}

theoretical_overlap_mi_genes <- function(d_current, d_reference, gamma, r = 10) {
  gamma_cur <- theoretical_gamma_genes_single(d_current, gamma, r = r)
  gamma_ref <- theoretical_gamma_genes_single(d_reference, gamma, r = r)
  gamma_cur[!is.finite(gamma_cur)] <- 0
  gamma_ref[!is.finite(gamma_ref)] <- 0
  gamma_double <- pmin(pmax(gamma_cur * gamma_ref, 0), 1 - 1e-12)
  list(
    mi = overlap_mi_from_gamma(gamma_double),
    gamma_current = gamma_cur,
    gamma_reference = gamma_ref,
    gamma_double = gamma_double
  )
}

theoretical_overlap_mi_from_spike_tables <- function(
  current_spikes,
  reference_spikes,
  side = c("cells", "genes"),
  gamma_current = NULL,
  gamma_reference = NULL,
  r = 10
) {
  side <- match.arg(side)
  if (!"spike_strength_d" %in% names(current_spikes)) {
    stop("current_spikes must contain a spike_strength_d column from fit_bbp_spikes().", call. = FALSE)
  }
  if (!"spike_strength_d" %in% names(reference_spikes)) {
    stop("reference_spikes must contain a spike_strength_d column from fit_bbp_spikes().", call. = FALSE)
  }

  d_current <- current_spikes$spike_strength_d
  d_reference <- reference_spikes$spike_strength_d
  if (!length(d_current) || !length(d_reference)) {
    gamma_zero <- rep(0, r)
    return(list(
      mi = 0,
      gamma_current = gamma_zero,
      gamma_reference = gamma_zero,
      gamma_double = gamma_zero
    ))
  }
  if (is.null(gamma_current)) {
    if (!"gamma" %in% names(current_spikes) || !length(current_spikes$gamma)) {
      stop("Provide gamma_current or include gamma in current_spikes.", call. = FALSE)
    }
    gamma_current <- current_spikes$gamma[1]
  }
  if (is.null(gamma_reference)) {
    if (!"gamma" %in% names(reference_spikes) || !length(reference_spikes$gamma)) {
      stop("Provide gamma_reference or include gamma in reference_spikes.", call. = FALSE)
    }
    gamma_reference <- reference_spikes$gamma[1]
  }

  if (side == "cells") {
    return(theoretical_overlap_mi_cells(d_current, d_reference, gamma_current, gamma_reference, r = r))
  }
  if (!isTRUE(all.equal(gamma_current, gamma_reference))) {
    warning("Gene-side formula currently uses gamma_current; gamma_reference is ignored.", call. = FALSE)
  }
  theoretical_overlap_mi_genes(d_current, d_reference, gamma_current, r = r)
}
