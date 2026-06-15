#!/usr/bin/env Rscript

source("scripts/core/mp_bbp.R")

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

eigenvalue_file <- arg_value("--eigenvalues")
out_dir <- arg_value("--out-dir", "outputs/bbp_fit")
dataset <- arg_value("--dataset", "dataset")
n_cells <- as.integer(arg_value("--n-cells", NA))
n_features <- as.integer(arg_value("--n-features", NA))
fit_type <- arg_value("--fit-type", "median_fixed_gamma")
edge_buffer <- as.numeric(arg_value("--edge-buffer", "1.01"))

if (is.null(eigenvalue_file) || !nzchar(eigenvalue_file)) {
  stop("Usage: run_bbp_fit_from_eigenvalues.R --eigenvalues=<csv> --out-dir=<dir> [--dataset=name] [--n-cells=N] [--n-features=P] [--fit-type=median_fixed_gamma|fixed_gamma|free_gamma]", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
eig_df <- read.csv(eigenvalue_file)
if (!"eigenvalue" %in% names(eig_df)) stop("CSV must contain an 'eigenvalue' column.", call. = FALSE)

if (is.na(n_cells)) {
  if ("n_cells" %in% names(eig_df)) {
    available_n <- eig_df$n_cells[is.finite(eig_df$n_cells)]
    if (length(available_n)) {
      n_cells <- max(available_n, na.rm = TRUE)
      eig_df <- eig_df[is.na(eig_df$n_cells) | eig_df$n_cells == n_cells, , drop = FALSE]
    }
  }
  if (is.na(n_cells)) {
    stop("Provide --n-cells or include finite n_cells in the CSV.", call. = FALSE)
  }
}

if (is.na(n_features)) {
  if ("n_features" %in% names(eig_df)) {
    available_p <- eig_df$n_features[is.finite(eig_df$n_features)]
    if (length(available_p)) {
      n_features <- max(available_p, na.rm = TRUE)
    }
  }
  if (is.na(n_features) && "n_genes" %in% names(eig_df)) {
    available_p <- eig_df$n_genes[is.finite(eig_df$n_genes)]
    if (length(available_p)) {
      n_features <- max(available_p, na.rm = TRUE)
    }
  }
  if (is.na(n_features)) {
    n_features <- sum(is.finite(eig_df$eigenvalue) & eig_df$eigenvalue > 0)
  }
}

if (is.na(n_features) || n_features <= 0) {
  stop("Could not infer n_features. Provide --n-features.", call. = FALSE)
}

fit <- fit_bbp_spikes(
  eig_df$eigenvalue,
  n_cells = n_cells,
  n_features = n_features,
  fit_type = fit_type,
  edge_buffer = edge_buffer
)

bulk_fit <- fit$bulk_fit
bulk_fit$dataset <- dataset
bulk_fit$eigenvalue_file <- eigenvalue_file
bulk_fit <- bulk_fit[, c("dataset", "eigenvalue_file", setdiff(names(bulk_fit), c("dataset", "eigenvalue_file")))]

curve <- mp_curve(bulk_fit$gamma, bulk_fit$sigma2)
spikes <- fit$spikes
if (nrow(spikes)) {
  spikes$dataset <- dataset
  spikes <- spikes[, c("dataset", setdiff(names(spikes), "dataset"))]
}

write.csv(bulk_fit, file.path(out_dir, "bbp_bulk_fit.csv"), row.names = FALSE)
write.csv(spikes, file.path(out_dir, "bbp_spikes.csv"), row.names = FALSE)
write.csv(curve, file.path(out_dir, "mp_bulk_curve.csv"), row.names = FALSE)

cat("Wrote BBP/MP fit outputs to", out_dir, "\n")
print(bulk_fit)
if (nrow(spikes)) print(utils::head(spikes, 20)) else cat("No BBP spikes above edge buffer.\n")
