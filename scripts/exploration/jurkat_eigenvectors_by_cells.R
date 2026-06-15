library(Seurat)
library(Matrix)
library(ggplot2)

set.seed(1)

data_dir <- "data/Jurkat/sample_filtered_feature_bc_matrix"
out_dir <- "outputs/exploration/jurkat_eigenvectors_by_cells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_variable_genes <- 2000
cell_counts_requested <- c(500, 1000, 2000, 5000, 10000, 25000, 50000, Inf)
top_k <- 50
subspace_sizes <- c(1, 2, 5, 10, 20, 50)
selected_ranks <- c(1, 2, 3, 4, 5, 10, 20, 50)

alignment_file <- file.path(out_dir, "jurkat_eigenvector_alignment_by_cells.csv")
subspace_file <- file.path(out_dir, "jurkat_eigenvector_subspace_alignment_by_cells.csv")
match_file <- file.path(out_dir, "jurkat_eigenvector_rank_match_by_cells.csv")
reference_file <- file.path(out_dir, "jurkat_full_reference_eigendecomp_top50.rds")
progress_file <- file.path(out_dir, "progress.log")

log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = progress_file, append = TRUE)
}

compute_cov_eigendecomp <- function(X, top_k) {
  n <- ncol(X)
  gene_means <- Matrix::rowMeans(X)
  cov_mat <- as.matrix(Matrix::tcrossprod(X))
  cov_mat <- (cov_mat - n * tcrossprod(gene_means)) / (n - 1)
  cov_mat <- (cov_mat + t(cov_mat)) / 2

  eig <- eigen(cov_mat, symmetric = TRUE)
  eig$values[eig$values < 0 & eig$values > -1e-8] <- 0

  keep <- seq_len(min(top_k, length(eig$values)))
  list(
    values = eig$values[keep],
    vectors = eig$vectors[, keep, drop = FALSE]
  )
}

label_cells <- function(x) {
  paste0("n = ", format(x, big.mark = ",", scientific = FALSE, trim = TRUE))
}

clamp_unit <- function(x) {
  pmin(pmax(x, 0), 1)
}

log_step("Loading 10x matrix from ", data_dir)
counts <- Read10X(data.dir = data_dir)
if (is.list(counts)) {
  log_step("Read10X returned feature types: ", paste(names(counts), collapse = ", "))
  if ("Gene Expression" %in% names(counts)) {
    counts <- counts[["Gene Expression"]]
  } else {
    counts <- counts[[1]]
    log_step("Gene Expression feature type not found; using first matrix")
  }
}

seu <- CreateSeuratObject(counts)
rm(counts)
gc()

log_step("Normalizing and selecting ", n_variable_genes, " HVGs")
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

n_genes <- nrow(X_full)
n_cells_total <- ncol(X_full)
cell_counts <- pmin(cell_counts_requested, n_cells_total)
cell_counts <- unique(as.integer(cell_counts))

log_step("Genes used: ", n_genes)
log_step("Total cells available: ", n_cells_total)
log_step("Cell counts analyzed: ", paste(cell_counts, collapse = ", "))
log_step("Top eigenvectors analyzed: ", top_k)

cell_order <- sample.int(n_cells_total)

if (file.exists(reference_file)) {
  log_step("Loading full reference eigendecomposition from ", reference_file)
  full_ref <- readRDS(reference_file)
} else {
  log_step("Computing full-cell reference eigendecomposition")
  full_ref <- compute_cov_eigendecomp(X_full, top_k)
  saveRDS(full_ref, reference_file)
  log_step("Saved full reference eigendecomposition to ", reference_file)
}

alignment_list <- vector("list", length(cell_counts))
subspace_list <- vector("list", length(cell_counts))
match_list <- vector("list", length(cell_counts))

for (i in seq_along(cell_counts)) {
  n_cells <- cell_counts[i]
  log_step("Computing eigenvectors for n_cells = ", n_cells)

  cell_idx <- if (n_cells == n_cells_total) {
    seq_len(n_cells_total)
  } else {
    sort(cell_order[seq_len(n_cells)])
  }

  eig <- if (n_cells == n_cells_total) {
    full_ref
  } else {
    compute_cov_eigendecomp(X_full[, cell_idx, drop = FALSE], top_k)
  }

  n_keep <- min(ncol(eig$vectors), ncol(full_ref$vectors))
  overlap <- abs(crossprod(eig$vectors[, seq_len(n_keep), drop = FALSE], full_ref$vectors[, seq_len(n_keep), drop = FALSE]))

  alignment_list[[i]] <- data.frame(
    n_cells = n_cells,
    rank = seq_len(n_keep),
    eigenvalue = eig$values[seq_len(n_keep)],
    full_eigenvalue = full_ref$values[seq_len(n_keep)],
    abs_cosine_to_full_same_rank = diag(overlap)
  )

  subspace_list[[i]] <- do.call(rbind, lapply(subspace_sizes[subspace_sizes <= n_keep], function(top_m) {
    singular_values <- svd(
      crossprod(
        eig$vectors[, seq_len(top_m), drop = FALSE],
        full_ref$vectors[, seq_len(top_m), drop = FALSE]
      ),
      nu = 0,
      nv = 0
    )$d

    data.frame(
      n_cells = n_cells,
      top_m = top_m,
      mean_canonical_correlation = mean(singular_values),
      min_canonical_correlation = min(singular_values)
    )
  }))

  match_list[[i]] <- do.call(rbind, lapply(seq_len(n_keep), function(rank) {
    best_ref_rank <- which.max(overlap[rank, ])
    data.frame(
      n_cells = n_cells,
      rank = rank,
      best_ref_rank = best_ref_rank,
      best_abs_cosine_to_full = overlap[rank, best_ref_rank],
      same_rank_abs_cosine_to_full = overlap[rank, rank]
    )
  }))

  gc()
}

alignment_df <- do.call(rbind, alignment_list)
subspace_df <- do.call(rbind, subspace_list)
match_df <- do.call(rbind, match_list)

alignment_df$abs_cosine_to_full_same_rank <- clamp_unit(alignment_df$abs_cosine_to_full_same_rank)
subspace_df$mean_canonical_correlation <- clamp_unit(subspace_df$mean_canonical_correlation)
subspace_df$min_canonical_correlation <- clamp_unit(subspace_df$min_canonical_correlation)
match_df$best_abs_cosine_to_full <- clamp_unit(match_df$best_abs_cosine_to_full)
match_df$same_rank_abs_cosine_to_full <- clamp_unit(match_df$same_rank_abs_cosine_to_full)

alignment_df$n_cells_label <- factor(label_cells(alignment_df$n_cells), levels = label_cells(cell_counts))
subspace_df$n_cells_label <- factor(label_cells(subspace_df$n_cells), levels = label_cells(cell_counts))
match_df$n_cells_label <- factor(label_cells(match_df$n_cells), levels = label_cells(cell_counts))

write.csv(alignment_df, alignment_file, row.names = FALSE)
write.csv(subspace_df, subspace_file, row.names = FALSE)
write.csv(match_df, match_file, row.names = FALSE)

theme_set(theme_classic(base_size = 13))

alignment_plot <- ggplot(alignment_df, aes(rank, abs_cosine_to_full_same_rank, color = n_cells_label)) +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Jurkat HVG eigenvector alignment to all-cell reference",
    subtitle = paste0("Same-rank absolute cosine similarity; top ", top_k, " covariance eigenvectors"),
    x = "Eigenvector rank",
    y = "Absolute cosine similarity",
    color = "Cells"
  )

selected_alignment <- subset(alignment_df, rank %in% selected_ranks)
selected_alignment$rank_label <- factor(paste0("PC", selected_alignment$rank), levels = paste0("PC", selected_ranks))

trajectory_plot <- ggplot(selected_alignment, aes(n_cells, abs_cosine_to_full_same_rank, color = rank_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Selected Jurkat eigenvector stability trajectories",
    x = "Number of cells",
    y = "Absolute cosine similarity to all-cell reference",
    color = "Rank"
  )

subspace_plot <- ggplot(subspace_df, aes(n_cells, mean_canonical_correlation, color = factor(top_m))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Jurkat top-m eigenspace alignment to all-cell reference",
    x = "Number of cells",
    y = "Mean canonical correlation",
    color = "Top m"
  )

rank_match_plot <- ggplot(match_df, aes(rank, best_ref_rank, fill = best_abs_cosine_to_full)) +
  geom_tile() +
  geom_abline(slope = 1, intercept = 0, color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(limits = c(0, 1), option = "magma") +
  facet_wrap(~ n_cells_label, ncol = 2) +
  coord_equal() +
  labs(
    title = "Best-matching all-cell reference rank",
    subtitle = "Diagonal structure indicates stable rank ordering",
    x = "Subset eigenvector rank",
    y = "Best reference rank",
    fill = "Abs cosine"
  )

ggsave(
  file.path(out_dir, "jurkat_eigenvector_alignment_by_cells.png"),
  alignment_plot,
  width = 9,
  height = 6,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_selected_eigenvector_trajectories.png"),
  trajectory_plot,
  width = 9,
  height = 6,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_subspace_alignment_by_cells.png"),
  subspace_plot,
  width = 9,
  height = 6,
  dpi = 220
)

ggsave(
  file.path(out_dir, "jurkat_eigenvector_best_reference_rank_by_cells.png"),
  rank_match_plot,
  width = 10,
  height = 12,
  dpi = 220
)

log_step("Wrote outputs to ", out_dir)
print(head(alignment_df, 20))
print(subspace_df)
