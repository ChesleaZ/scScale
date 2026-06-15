#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
pkg_dir <- normalizePath(file.path(script_dir, "..", ".."))
for (f in sort(list.files(file.path(pkg_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

cur_file <- arg_value("--cur-eigenvalues")
ref_file <- arg_value("--ref-eigenvalues", cur_file)
cur_counts_file <- arg_value("--cur-counts")
ref_counts_file <- arg_value("--ref-counts", cur_counts_file)
out_dir <- arg_value("--out-dir", "outputs/scscale_eigenvalues")
side <- arg_value("--side", "cells")
n_cur <- as.integer(arg_value("--n-cur"))
p_cur <- as.integer(arg_value("--p-cur"))
n_ref <- as.integer(arg_value("--n-ref", as.character(n_cur)))
p_ref <- as.integer(arg_value("--p-ref", as.character(p_cur)))
r <- as.integer(arg_value("--r", "10"))
R <- as.integer(arg_value("--mp-R", "50"))
p_sim <- as.integer(arg_value("--mp-p-sim", "800"))
seed <- as.integer(arg_value("--seed", "1"))
n_features <- as.integer(arg_value("--n-features", "2000"))
min_cells <- as.integer(arg_value("--min-cells", "10"))
transform <- arg_value("--transform", "pearson")

if (is.null(cur_file) && is.null(cur_counts_file)) {
  stop(
    paste(
      "Usage:",
      "Rscript packages/scScale/inst/scripts/run_from_eigenvalues.R",
      "--cur-eigenvalues=<csv> [--ref-eigenvalues=<csv>]",
      "--n-cur=<N> --p-cur=<P> [--n-ref=<N> --p-ref=<P>]",
      "or --cur-counts=<csv> [--ref-counts=<csv>]",
      "[--side=cells|genes] [--r=10] [--out-dir=<dir>]"
    ),
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(cur_file)) {
  if (is.na(n_cur) || is.na(p_cur)) {
    stop("--n-cur and --p-cur are required when using precomputed eigenvalue CSVs.", call. = FALSE)
  }
  cur_ev <- find_eigenvalues(read_eigenvalues_csv(cur_file), n = n_cur, p = p_cur)
  ref_ev <- find_eigenvalues(read_eigenvalues_csv(ref_file), n = n_ref, p = p_ref)
} else {
  cur_ev <- find_eigenvalues(
    read_counts_csv(cur_counts_file),
    n_features = n_features,
    min_cells = min_cells,
    transform = transform
  )
  ref_ev <- find_eigenvalues(
    read_counts_csv(ref_counts_file),
    n_features = n_features,
    min_cells = min_cells,
    transform = transform
  )
}

fit <- mi_theory_from_eigenvalues(
  cur_eigenvalues = cur_ev,
  ref_eigenvalues = ref_ev,
  side = side,
  r = r,
  R = R,
  p_sim = p_sim,
  seed = seed
)

summary_df <- data.frame(
  cur_eigenvalue_file = cur_file %||% NA_character_,
  ref_eigenvalue_file = ref_file %||% NA_character_,
  cur_counts_file = cur_counts_file %||% NA_character_,
  ref_counts_file = ref_counts_file %||% NA_character_,
  side = side,
  r = r,
  n_cur = fit$cur$n,
  p_cur = fit$cur$p,
  c_cur = fit$cur$c_ratio,
  n_ref = fit$ref$n,
  p_ref = fit$ref$p,
  c_ref = fit$ref$c_ratio,
  noise_method = fit$noise$method,
  noise_scale_ref = fit$noise$scale,
  mp_upper_cur = fit$cur$mp_upper,
  mp_upper_ref = fit$ref$mp_upper,
  n_spikes_cur = fit$cur$n_spikes,
  n_spikes_ref = fit$ref$n_spikes,
  mi_theory = fit$mi
)

gamma_df <- data.frame(
  rank = seq_len(r),
  gamma_cur = fit$gamma_cur,
  gamma_ref = fit$gamma_ref,
  gamma_double = fit$gamma_double
)

utils::write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(gamma_df, file.path(out_dir, "theoretical_gamma.csv"), row.names = FALSE)
utils::write.csv(fit$cur$spike_table, file.path(out_dir, "current_spikes.csv"), row.names = FALSE)
utils::write.csv(fit$ref$spike_table, file.path(out_dir, "reference_spikes.csv"), row.names = FALSE)

cat("Wrote scScale eigenvalue outputs to", out_dir, "\n")
print(summary_df)
