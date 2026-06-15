library(Seurat)
library(Matrix)
library(ggplot2)

set.seed(1)

data_dir <- "data/Jurkat/sample_filtered_feature_bc_matrix"
out_dir <- "outputs/exploration/jurkat_eigenvalues_by_cells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_variable_genes <- 2000
cell_counts_requested <- c(500, 1000, 2000, 5000, 10000, 25000, 50000, Inf)
clip_at <- 2

compute_cov_eigenvalues <- function(X) {
  n <- ncol(X)
  gene_means <- Matrix::rowMeans(X)
  cov_mat <- as.matrix(Matrix::tcrossprod(X))
  cov_mat <- (cov_mat - n * tcrossprod(gene_means)) / (n - 1)
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  eigvals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eigvals[eigvals < 0 & eigvals > -1e-8] <- 0
  sort(eigvals, decreasing = TRUE)
}

message("Loading 10x matrix from ", data_dir)
counts <- Read10X(data.dir = data_dir)
seu <- CreateSeuratObject(counts)

message("Normalizing and selecting variable genes")
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = n_variable_genes,
  verbose = FALSE
)

genes_use <- VariableFeatures(seu)
X_full <- GetAssayData(seu, assay = DefaultAssay(seu), layer = "data")[genes_use, ]

n_cells_total <- ncol(X_full)
cell_counts <- pmin(cell_counts_requested, n_cells_total)
cell_counts <- unique(as.integer(cell_counts))

message("Genes used: ", nrow(X_full))
message("Total cells available: ", n_cells_total)
message("Cell counts analyzed: ", paste(cell_counts, collapse = ", "))

eig_list <- vector("list", length(cell_counts))
summary_list <- vector("list", length(cell_counts))

for (i in seq_along(cell_counts)) {
  n_cells <- cell_counts[i]
  message("Computing eigenvalues for n_cells = ", n_cells)

  cell_idx <- if (n_cells == n_cells_total) {
    seq_len(n_cells_total)
  } else {
    sort(sample.int(n_cells_total, n_cells))
  }

  eigvals <- compute_cov_eigenvalues(X_full[, cell_idx, drop = FALSE])

  eig_list[[i]] <- data.frame(
    n_cells = n_cells,
    rank = seq_along(eigvals),
    eigenvalue = eigvals,
    eigenvalue_clipped = pmin(eigvals, clip_at),
    clipped = eigvals > clip_at
  )

  summary_list[[i]] <- data.frame(
    n_cells = n_cells,
    n_genes = nrow(X_full),
    max_eigenvalue = max(eigvals),
    median_eigenvalue = median(eigvals),
    mean_eigenvalue = mean(eigvals),
    n_eigenvalues_gt_2 = sum(eigvals > clip_at)
  )
}

eig_df <- do.call(rbind, eig_list)
summary_df <- do.call(rbind, summary_list)

eig_df$n_cells_label <- factor(
  eig_df$n_cells,
  levels = cell_counts,
  labels = paste0("n = ", format(cell_counts, big.mark = ","))
)

write.csv(
  eig_df,
  file.path(out_dir, "jurkat_eigenvalues_by_cells.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(out_dir, "jurkat_eigenvalues_by_cells_summary.csv"),
  row.names = FALSE
)

theme_set(theme_classic(base_size = 13))

hist_plot <- ggplot(eig_df, aes(eigenvalue_clipped)) +
  geom_histogram(bins = 70, fill = "grey70", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = clip_at, color = "#C43C39", linetype = "dashed") +
  facet_wrap(~ n_cells_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Jurkat eigenvalue distributions by number of cells",
    subtitle = paste0("Eigenvalues larger than ", clip_at, " are clipped to ", clip_at, " for plotting"),
    x = "Covariance eigenvalue, clipped",
    y = "Count"
  )

density_plot <- ggplot(eig_df, aes(eigenvalue_clipped, color = n_cells_label)) +
  geom_density(linewidth = 0.9, adjust = 0.8) +
  geom_vline(xintercept = clip_at, color = "#C43C39", linetype = "dashed") +
  labs(
    title = "Overlay of clipped Jurkat eigenvalue distributions",
    subtitle = paste0("Eigenvalues larger than ", clip_at, " are clipped to ", clip_at),
    x = "Covariance eigenvalue, clipped",
    y = "Density",
    color = "Cells"
  )

tail_plot <- ggplot(summary_df, aes(n_cells, n_eigenvalues_gt_2)) +
  geom_line(linewidth = 0.9, color = "grey25") +
  geom_point(size = 2.2, color = "#1F78B4") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Number of eigenvalues greater than 2",
    x = "Number of cells",
    y = "Eigenvalues > 2"
  )

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_distributions_by_cells_clipped.png"),
  hist_plot,
  width = 11,
  height = 12,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_eigenvalue_distributions_by_cells_overlay.png"),
  density_plot,
  width = 9,
  height = 6,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_eigenvalues_gt_2_by_cells.png"),
  tail_plot,
  width = 8,
  height = 5,
  dpi = 220
)

print(summary_df)
message("Wrote outputs to ", out_dir)
