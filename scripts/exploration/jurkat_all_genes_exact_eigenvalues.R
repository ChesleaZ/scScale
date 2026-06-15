library(Matrix)
library(ggplot2)

set.seed(1)

data_dir <- "data/Jurkat/sample_filtered_feature_bc_matrix"
out_dir <- "outputs/exploration/jurkat_all_genes_exact"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cov_file <- file.path(out_dir, "jurkat_all_genes_covariance.rds")
eig_rds_file <- file.path(out_dir, "jurkat_all_genes_eigenvalues.rds")
eig_csv_file <- file.path(out_dir, "jurkat_all_genes_eigenvalues.csv")
summary_file <- file.path(out_dir, "jurkat_all_genes_summary.csv")
progress_file <- file.path(out_dir, "progress.log")

log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = progress_file, append = TRUE)
}

read_10x_sparse <- function(data_dir) {
  matrix_file <- file.path(data_dir, "matrix.mtx.gz")
  features_file <- file.path(data_dir, "features.tsv.gz")
  barcodes_file <- file.path(data_dir, "barcodes.tsv.gz")

  counts <- Matrix::readMM(gzfile(matrix_file))
  counts <- as(counts, "dgCMatrix")

  features <- read.delim(
    gzfile(features_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  barcodes <- read.delim(
    gzfile(barcodes_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(counts) <- make.unique(features[[2]])
  colnames(counts) <- barcodes[[1]]
  counts
}

normalize_log_sparse <- function(counts, scale_factor = 1e4) {
  cell_sums <- Matrix::colSums(counts)
  if (any(cell_sums == 0)) {
    stop("Found zero-count cells; cannot LogNormalize by cell depth.")
  }

  x <- counts
  x@x <- log1p(x@x / rep.int(cell_sums, diff(x@p)) * scale_factor)
  x
}

compute_covariance <- function(x) {
  n <- ncol(x)
  gene_means <- Matrix::rowMeans(x)

  log_step("Computing sparse tcrossprod for ", nrow(x), " genes x ", n, " cells")
  cov_mat <- as.matrix(Matrix::tcrossprod(x))

  log_step("Centering covariance matrix")
  cov_mat <- (cov_mat - n * tcrossprod(gene_means)) / (n - 1)
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  cov_mat
}

if (file.exists(eig_csv_file) && file.exists(eig_rds_file)) {
  log_step("Eigenvalues already exist; skipping computation")
  eigvals <- readRDS(eig_rds_file)
} else {
  if (file.exists(cov_file)) {
    log_step("Loading checkpointed covariance from ", cov_file)
    cov_mat <- readRDS(cov_file)
  } else {
    log_step("Reading full 10x count matrix")
    counts <- read_10x_sparse(data_dir)
    log_step(
      "Loaded counts: ",
      nrow(counts), " genes x ", ncol(counts), " cells; nnz = ", length(counts@x)
    )

    log_step("Log-normalizing sparse count matrix")
    x <- normalize_log_sparse(counts)
    rm(counts)
    gc()

    cov_mat <- compute_covariance(x)
    rm(x)
    gc()

    log_step("Saving covariance checkpoint to ", cov_file)
    saveRDS(cov_mat, cov_file, compress = FALSE)
    gc()
  }

  log_step("Computing exact eigenvalues of dense covariance")
  eigvals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eigvals[eigvals < 0 & eigvals > -1e-8] <- 0
  eigvals <- sort(eigvals, decreasing = TRUE)

  log_step("Saving eigenvalue checkpoint to ", eig_rds_file)
  saveRDS(eigvals, eig_rds_file, compress = FALSE)

  eig_df <- data.frame(
    rank = seq_along(eigvals),
    eigenvalue = eigvals,
    log10_eigenvalue = log10(pmax(eigvals, .Machine$double.eps))
  )
  write.csv(eig_df, eig_csv_file, row.names = FALSE)

  rm(cov_mat)
  gc()
}

summary_df <- data.frame(
  metric = c(
    "n_genes",
    "min_eigenvalue",
    "median_eigenvalue",
    "mean_eigenvalue",
    "max_eigenvalue",
    "n_eigenvalues_gt_0.5",
    "n_eigenvalues_gt_2"
  ),
  value = c(
    length(eigvals),
    min(eigvals),
    median(eigvals),
    mean(eigvals),
    max(eigvals),
    sum(eigvals > 0.5),
    sum(eigvals > 2)
  )
)
write.csv(summary_df, summary_file, row.names = FALSE)

eig_df <- data.frame(
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  log10_eigenvalue = log10(pmax(eigvals, .Machine$double.eps)),
  eigenvalue_clip_0_5 = pmin(eigvals, 0.5),
  eigenvalue_clip_2 = pmin(eigvals, 2)
)

theme_set(theme_classic(base_size = 14))

p_clip_0_5 <- ggplot(eig_df, aes(eigenvalue_clip_0_5)) +
  geom_histogram(bins = 220, fill = "grey72", color = "white", linewidth = 0.12) +
  geom_vline(xintercept = 0.5, color = "#C43C39", linetype = "dashed") +
  labs(
    title = "All-gene Jurkat eigenvalue distribution",
    subtitle = "All cells, all genes; eigenvalues clipped at 0.5",
    x = "Covariance eigenvalue, clipped",
    y = "Count"
  )

p_clip_2 <- ggplot(eig_df, aes(eigenvalue_clip_2)) +
  geom_histogram(bins = 220, fill = "grey72", color = "white", linewidth = 0.12) +
  geom_vline(xintercept = 2, color = "#C43C39", linetype = "dashed") +
  labs(
    title = "All-gene Jurkat eigenvalue distribution",
    subtitle = "All cells, all genes; eigenvalues clipped at 2",
    x = "Covariance eigenvalue, clipped",
    y = "Count"
  )

p_scree <- ggplot(eig_df, aes(rank, eigenvalue)) +
  geom_line(linewidth = 0.4, color = "grey35") +
  scale_y_log10() +
  labs(
    title = "All-gene Jurkat covariance spectrum",
    x = "Eigenvalue rank",
    y = "Eigenvalue"
  )

ggsave(
  file.path(out_dir, "jurkat_all_genes_eigenvalue_distribution_clip_0_5.png"),
  p_clip_0_5,
  width = 8,
  height = 5,
  dpi = 220
)
ggsave(
  file.path(out_dir, "jurkat_all_genes_eigenvalue_distribution_clip_2.png"),
  p_clip_2,
  width = 8,
  height = 5,
  dpi = 220
)
ggsave(
  file.path(out_dir, "jurkat_all_genes_log_scree.png"),
  p_scree,
  width = 8,
  height = 5,
  dpi = 220
)

print(summary_df)
log_step("Done")
