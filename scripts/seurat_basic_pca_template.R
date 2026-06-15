#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
dataset_name <- if (length(args) >= 1) args[[1]] else "jurkat"
input_type <- if (length(args) >= 2) args[[2]] else "tenx_mtx"
input_path <- if (length(args) >= 3) args[[3]] else "data/Jurkat/sample_filtered_feature_bc_matrix"
out_dir <- if (length(args) >= 4) args[[4]] else file.path("outputs", "processed", dataset_name)

normalization_method <- "LogNormalize"
scale_factor <- 10000
hvg_method <- "vst"
n_hvg <- 2000
n_pcs_requested <- if (length(args) >= 5) as.integer(args[[5]]) else 50
umap_dims_requested <- 30

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("scripts/core/mp_bbp.R")

message("Dataset: ", dataset_name)
message("Input: ", input_type, " ", input_path)
message("Output: ", out_dir)

read_counts <- function(input_type, input_path) {
  if (input_type == "tenx_mtx") {
    counts <- Read10X(data.dir = input_path)
    if (is.list(counts)) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }
    return(counts)
  }
  stop("Unsupported input_type for this template: ", input_type, call. = FALSE)
}

counts <- read_counts(input_type, input_path)

obj <- CreateSeuratObject(counts = counts, project = dataset_name)
obj <- NormalizeData(
  obj,
  normalization.method = normalization_method,
  scale.factor = scale_factor,
  verbose = FALSE
)
obj <- FindVariableFeatures(
  obj,
  selection.method = hvg_method,
  nfeatures = n_hvg,
  verbose = FALSE
)

hvg <- VariableFeatures(obj)
n_pcs <- min(n_pcs_requested, length(hvg) - 1L, ncol(obj) - 1L)
umap_dims <- seq_len(min(umap_dims_requested, n_pcs))

obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = n_pcs, verbose = FALSE)

eigvals <- obj[["pca"]]@stdev^2
eig_df <- data.frame(
  dataset = dataset_name,
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  log10_eigenvalue = log10(pmax(eigvals, .Machine$double.eps))
)
write.csv(eig_df, file.path(out_dir, "pca_eigenvalues.csv"), row.names = FALSE)

mp_fit <- fit_mp_median(
  eigvals,
  n_cells = ncol(obj),
  n_features = length(hvg),
  edge_buffer = 1.01
)
write.csv(mp_fit, file.path(out_dir, "mp_fit_parameters.csv"), row.names = FALSE)

curve_df <- mp_curve(mp_fit$gamma, mp_fit$sigma2, n = 1000, positive_only = TRUE)
write.csv(curve_df, file.path(out_dir, "mp_fit_curve.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 12))

scree_plot <- ggplot(eig_df, aes(rank, eigenvalue)) +
  geom_line(color = "grey35", linewidth = 0.45) +
  geom_point(size = 0.9, color = "grey20") +
  geom_hline(yintercept = mp_fit$mp_upper, color = "#1F78B4", linetype = "dashed") +
  scale_y_log10() +
  labs(
    title = paste0(dataset_name, ": PCA eigenvalues"),
    subtitle = "Dashed blue: median fixed-gamma MP upper edge",
    x = "PC rank",
    y = "Eigenvalue"
  )

bulk_df <- eig_df[eig_df$eigenvalue > 1e-10 & eig_df$eigenvalue <= mp_fit$mp_upper * 1.25, ]
mp_plot <- ggplot(bulk_df, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey78", color = "white", linewidth = 0.15) +
  geom_line(data = curve_df, aes(eigenvalue, density), color = "#1F78B4", linewidth = 0.9) +
  geom_vline(xintercept = mp_fit$mp_upper, color = "#C43C39", linetype = "dashed") +
  labs(
    title = paste0(dataset_name, ": residual PCA bulk vs MP fit"),
    subtitle = "Median fixed-gamma MP fit",
    x = "PCA eigenvalue",
    y = "Density"
  )

ggsave(file.path(out_dir, "pca_eigenvalue_scree.png"), scree_plot, width = 7, height = 5, dpi = 220)
ggsave(file.path(out_dir, "pca_eigenvalue_mp_fit.png"), mp_plot, width = 7, height = 5, dpi = 220)

obj <- RunUMAP(obj, dims = umap_dims, verbose = FALSE)
umap_plot <- DimPlot(obj, reduction = "umap") +
  ggtitle(paste0(dataset_name, ": default UMAP"))
ggsave(file.path(out_dir, "umap_default.png"), umap_plot, width = 7, height = 5, dpi = 220)

metadata <- list(
  dataset = dataset_name,
  input_type = input_type,
  input_path = input_path,
  n_cells = ncol(obj),
  n_features_total = nrow(obj),
  normalization_method = normalization_method,
  scale_factor = scale_factor,
  hvg_method = hvg_method,
  n_hvg = length(hvg),
  n_pcs = n_pcs,
  umap_dims = umap_dims,
  mp_fit = mp_fit
)

obj_to_save <- DietSeurat(
  obj,
  assays = DefaultAssay(obj),
  features = hvg,
  layers = c("counts", "data"),
  dimreducs = c("pca", "umap"),
  graphs = NULL,
  misc = TRUE
)
obj_to_save@misc$processing_metadata <- metadata

cache <- list(
  metadata = metadata,
  variable_features = hvg,
  pca_embeddings = Embeddings(obj, "pca"),
  pca_loadings = Loadings(obj, "pca"),
  pca_stdev = obj[["pca"]]@stdev,
  pca_eigenvalues = eig_df,
  umap_embeddings = Embeddings(obj, "umap"),
  mp_fit = mp_fit,
  mp_curve = curve_df
)

saveRDS(cache, file.path(out_dir, "pca_cache.rds"), compress = FALSE)
saveRDS(obj_to_save, file.path(out_dir, "seurat_basic_processed.rds"), compress = TRUE)
saveRDS(metadata, file.path(out_dir, "processing_metadata.rds"))
writeLines(hvg, file.path(out_dir, "variable_features.txt"))

message("Wrote reusable Seurat object and diagnostics to ", out_dir)
print(mp_fit)
