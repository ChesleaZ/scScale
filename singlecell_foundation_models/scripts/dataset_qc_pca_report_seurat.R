#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1]]
}

has_flag <- function(flag) flag %in% args

required_arg <- function(flag) {
  value <- get_arg(flag)
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument: ", flag, call. = FALSE)
  }
  value
}

dataset_name <- required_arg("--dataset-name")
input_type <- required_arg("--input-type")
input_path <- required_arg("--input-path")
out_dir <- required_arg("--out-dir")
technology <- get_arg("--technology", "unknown")
assay_label <- get_arg("--assay", "scRNA-seq")
min_genes <- as.integer(get_arg("--min-genes", "200"))
min_cells <- as.integer(get_arg("--min-cells", "3"))
target_sum <- as.numeric(get_arg("--target-sum", "10000"))
n_hvg <- as.integer(get_arg("--n-hvg", "2000"))
embedding_pcs_requested <- as.integer(get_arg("--embedding-pcs", "30"))
neighbors <- as.integer(get_arg("--neighbors", "15"))
cluster_resolution <- as.numeric(get_arg("--cluster-resolution", "1.0"))
max_cells_arg <- get_arg("--max-cells", NA)
max_cells <- if (is.na(max_cells_arg)) NA_integer_ else as.integer(max_cells_arg)
seed <- as.integer(get_arg("--seed", "0"))
transpose_table <- has_flag("--transpose-table")
save_rds <- has_flag("--save-rds")
python_bin <- get_arg("--python", "python")
stage_root <- get_arg("--stage-root", file.path(out_dir, "_staged_input"))
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]))) else getwd()
core_candidates <- c(
  file.path(getwd(), "scripts/core/mp_bbp.R"),
  file.path(dirname(dirname(script_dir)), "scripts/core/mp_bbp.R"),
  file.path(dirname(dirname(dirname(script_dir))), "scripts/core/mp_bbp.R")
)
core_mp <- core_candidates[file.exists(core_candidates)][1]
if (is.na(core_mp)) stop("Could not find scripts/core/mp_bbp.R", call. = FALSE)
source(core_mp)

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

fmt_num <- function(x, digits = 2) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  ax <- abs(as.numeric(x))
  if (ax >= 1e9) return(sprintf(paste0("%.", digits, "fB"), x / 1e9))
  if (ax >= 1e6) return(sprintf(paste0("%.", digits, "fM"), x / 1e6))
  if (ax >= 1e3) return(format(round(x), big.mark = ",", scientific = FALSE))
  if (abs(x - round(x)) < .Machine$double.eps^0.5) return(as.character(round(x)))
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}

choose_cells <- function(n, max_cells, seed) {
  if (is.na(max_cells) || n <= max_cells) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, max_cells))
}

read_counts <- function(input_type, input_path) {
  if (input_type == "10x_mtx") {
    counts <- Read10X(data.dir = input_path)
    if (is.list(counts)) {
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }
    return(counts)
  }
  if (input_type == "10x_h5") {
    counts <- Read10X_h5(filename = input_path)
    if (is.list(counts)) {
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }
    return(counts)
  }
  if (input_type %in% c("csv", "tsv")) {
    sep <- if (input_type == "csv") "," else "\t"
    mat <- as.matrix(read.table(input_path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE))
    if (transpose_table) mat <- t(mat)
    return(Matrix(t(mat), sparse = TRUE))
  }
  if (input_type %in% c("h5ad", "parse_zip")) {
    helper <- file.path(script_dir, "export_counts_for_seurat.py")
    if (!file.exists(helper)) {
      helper <- file.path(getwd(), "singlecell_foundation_models/scripts/export_counts_for_seurat.py")
    }
    dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
    cmd <- c(
      helper,
      "--input-type", input_type,
      "--input-path", input_path,
      "--out-dir", stage_root,
      "--seed", as.character(seed)
    )
    if (!is.na(max_cells)) cmd <- c(cmd, "--max-cells", as.character(max_cells))
    status <- system2(python_bin, cmd)
    if (!identical(status, 0L)) stop("Python staging failed for ", input_type, call. = FALSE)
    counts <- Read10X(data.dir = stage_root)
    if (is.list(counts)) {
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }
    return(counts)
  }
  stop("Unsupported input type: ", input_type, call. = FALSE)
}

get_layer <- function(obj, layer) {
  assay_name <- DefaultAssay(obj)
  out <- tryCatch(
    LayerData(obj, assay = assay_name, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay_name, slot = layer)
  )
  out
}

row_variance <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    means <- Matrix::rowMeans(x)
    means_sq <- Matrix::rowMeans(x ^ 2)
    vars <- means_sq - means ^ 2
  } else {
    vars <- apply(x, 1, var)
  }
  vars[!is.finite(vars)] <- -Inf
  vars
}

save_plot <- function(filename, plot, width = 7, height = 4.8) {
  ggsave(file.path(out_dir, filename), plot, width = width, height = height, dpi = 190, bg = "white")
  filename
}

table_html <- function(rows) {
  body <- paste(
    apply(rows, 1, function(z) sprintf("<tr><th>%s</th><td>%s</td></tr>", html_escape(z[[1]]), z[[2]])),
    collapse = "\n"
  )
  paste0("<table>", body, "</table>")
}

message("[", Sys.time(), "] Reading input")
counts <- read_counts(input_type, input_path)
if (!inherits(counts, "sparseMatrix")) counts <- Matrix(counts, sparse = TRUE)
counts <- as(counts, "dgCMatrix")

if (!is.na(max_cells) && ncol(counts) > max_cells && !(input_type %in% c("h5ad", "parse_zip"))) {
  keep <- choose_cells(ncol(counts), max_cells, seed)
  counts <- counts[, keep, drop = FALSE]
}

initial_cells <- ncol(counts)
initial_features <- nrow(counts)
umi <- Matrix::colSums(counts)
detected <- Matrix::colSums(counts > 0)

message("[", Sys.time(), "] Creating Seurat object")
obj <- CreateSeuratObject(counts = counts, project = dataset_name, min.cells = 0, min.features = 0)
obj[["nCount_report"]] <- umi[colnames(obj)]
obj[["nFeature_report"]] <- detected[colnames(obj)]

if (min_genes > 0) {
  obj <- subset(obj, subset = nFeature_RNA >= min_genes)
}
if (min_cells > 0) {
  counts_filtered <- get_layer(obj, "counts")
  keep_features <- Matrix::rowSums(counts_filtered > 0) >= min_cells
  obj <- subset(obj, features = rownames(counts_filtered)[keep_features])
}
filtered_cells <- ncol(obj)
filtered_features <- nrow(obj)

message("[", Sys.time(), "] NormalizeData")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = target_sum, verbose = FALSE)

message("[", Sys.time(), "] Naive top-variance HVG selection")
norm_data <- get_layer(obj, "data")
vars <- row_variance(norm_data)
top_n <- min(n_hvg, length(vars))
hvg <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
VariableFeatures(obj) <- hvg

message("[", Sys.time(), "] ScaleData and full PCA")
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
full_pcs <- min(length(hvg), ncol(obj) - 1L)
if (full_pcs < 2L) stop("Not enough cells/features after filtering for PCA", call. = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = full_pcs, approx = FALSE, verbose = FALSE)
eigvals <- obj[["pca"]]@stdev ^ 2
eig_df <- data.frame(
  rank = seq_along(eigvals),
  eigenvalue = eigvals,
  log10_eigenvalue = log10(pmax(eigvals, .Machine$double.xmin))
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
write.csv(
  data.frame(gene = hvg, variance = vars[hvg], rank = seq_along(hvg)),
  file.path(out_dir, "naive_hvg_top_variance.csv"),
  row.names = FALSE
)

theme_set(theme_classic(base_size = 12))
rank <- seq_along(umi)
rank_df <- data.frame(rank = rank, umi = sort(pmax(as.numeric(umi), 1), decreasing = TRUE))
rank_plot <- ggplot(rank_df, aes(rank, umi)) +
  geom_line(color = "#2454a6", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = paste0(dataset_name, ": UMI rank plot"), x = "Cell barcode rank", y = "UMI/counts per cell")
rank_file <- save_plot("umi_rank_plot.png", rank_plot)

umi_df <- data.frame(log10_umi = log10(pmax(as.numeric(umi), 1)))
umi_plot <- ggplot(umi_df, aes(log10_umi)) +
  geom_histogram(bins = 80, fill = "#18745a", color = "white", linewidth = 0.15) +
  labs(title = paste0(dataset_name, ": UMI distribution"), x = "log10(UMI/counts per cell)", y = "Cells")
umi_file <- save_plot("umi_distribution.png", umi_plot)

scree_plot <- ggplot(eig_df, aes(rank, eigenvalue)) +
  geom_line(color = "#17202a", linewidth = 0.45) +
  geom_point(size = 0.8, color = "#17202a") +
  geom_hline(yintercept = mp_fit$mp_upper, color = "#2454a6", linetype = "dashed") +
  scale_y_log10() +
  labs(title = paste0(dataset_name, ": PCA eigenvalues"), x = "PC rank", y = "Eigenvalue")
scree_file <- save_plot("pca_eigenvalue_scree.png", scree_plot)

bulk <- eig_df[eig_df$eigenvalue > 1e-10 & eig_df$eigenvalue <= max(mp_fit$mp_upper * 1.25, mp_fit$mp_upper), ]
mp_plot <- ggplot(bulk, aes(eigenvalue)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d9e1e8", color = "white", linewidth = 0.15) +
  geom_line(data = curve_df, aes(eigenvalue, density), color = "#2454a6", linewidth = 0.9) +
  geom_vline(xintercept = mp_fit$mp_upper, color = "#9f2a3b", linetype = "dashed") +
  labs(title = paste0(dataset_name, ": MP fit to residual eigenvalue bulk"), x = "PCA eigenvalue", y = "Density")
mp_file <- save_plot("pca_mp_bulk_fit.png", mp_plot)

message("[", Sys.time(), "] Neighbors, clusters, UMAP")
embedding_pcs <- min(embedding_pcs_requested, full_pcs)
dims_use <- seq_len(embedding_pcs)
obj <- FindNeighbors(obj, dims = dims_use, k.param = neighbors, verbose = FALSE)
obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
umap_plot <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = FALSE) +
  ggtitle(paste0(dataset_name, ": UMAP by Seurat cluster"))
umap_file <- save_plot("umap.png", umap_plot, width = 8, height = 6)

metadata <- list(
  dataset_name = dataset_name,
  technology = technology,
  assay = assay_label,
  input_type = input_type,
  input_path = input_path,
  cells_before_filtering = initial_cells,
  features_before_filtering = initial_features,
  cells_after_filtering = filtered_cells,
  features_after_filtering = filtered_features,
  total_umis = sum(as.numeric(umi)),
  median_umi_per_cell = median(as.numeric(umi)),
  mean_umi_per_cell = mean(as.numeric(umi)),
  median_features_per_cell = median(as.numeric(detected)),
  normalization_method = "Seurat::NormalizeData(LogNormalize)",
  scale_factor = target_sum,
  hvg_method = "naive_top_variance_on_log_normalized_data",
  n_hvg = length(hvg),
  full_pca_components = full_pcs,
  embedding_pcs = embedding_pcs,
  cluster_resolution = cluster_resolution,
  n_clusters = length(unique(obj$seurat_clusters)),
  mp_fit = mp_fit
)
write_json(metadata, file.path(out_dir, "report_metrics.json"), pretty = TRUE, auto_unbox = TRUE)

if (save_rds) {
  saveRDS(
    list(
      metadata = metadata,
      variable_features = hvg,
      pca_embeddings = Embeddings(obj, "pca"),
      pca_loadings = Loadings(obj, "pca"),
      pca_stdev = obj[["pca"]]@stdev,
      umap_embeddings = Embeddings(obj, "umap"),
      clusters = obj$seurat_clusters
    ),
    file.path(out_dir, "seurat_pca_cache.rds"),
    compress = TRUE
  )
  saveRDS(obj, file.path(out_dir, "seurat_processed.rds"), compress = TRUE)
}

qc_rows <- rbind(
  c("Dataset", html_escape(dataset_name)),
  c("Technology", html_escape(technology)),
  c("Assay", html_escape(assay_label)),
  c("Input", paste0("<code>", html_escape(input_path), "</code>")),
  c("Cells before filtering", fmt_num(initial_cells)),
  c("Features before filtering", fmt_num(initial_features)),
  c("Cells after filtering", fmt_num(filtered_cells)),
  c("Features after filtering", fmt_num(filtered_features)),
  c("Total UMIs/counts", fmt_num(sum(as.numeric(umi)))),
  c("Median UMI/cell", fmt_num(median(as.numeric(umi)))),
  c("Mean UMI/cell", fmt_num(mean(as.numeric(umi)))),
  c("Median detected features/cell", fmt_num(median(as.numeric(detected)))),
  c("Normalization", "Seurat LogNormalize"),
  c("HVG selection", "naive top variance on log-normalized data"),
  c("Highly variable genes", fmt_num(length(hvg))),
  c("Full PCA components stored", fmt_num(full_pcs)),
  c("PCs used for neighbors/UMAP", fmt_num(embedding_pcs)),
  c("Seurat clusters", fmt_num(length(unique(obj$seurat_clusters)))),
  c("MP fit", paste0(html_escape(mp_fit$fit_type), ", iterative quantile fit")),
  c("MP upper edge", fmt_num(mp_fit$mp_upper, digits = 4)),
  c("Eigenvalue outliers above MP edge", fmt_num(mp_fit$n_outliers_above_mp)),
  c("Top eigenvalue / MP edge", fmt_num(mp_fit$lambda1_over_mp_upper, digits = 2))
)

html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>%s Seurat QC + PCA Report</title>
  <style>
    body { margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17202a; background: #f7f9fb; }
    header { padding: 28px 34px; background: #fff; border-bottom: 1px solid #d9e1e8; }
    main { max-width: 1220px; margin: 0 auto; padding: 24px; }
    h1 { margin: 0 0 6px; font-size: 34px; }
    h2 { margin: 0 0 12px; font-size: 22px; }
    p { color: #667085; margin: 0; }
    code { background: #eef3f7; border: 1px solid #dbe5ec; padding: 1px 5px; border-radius: 5px; }
    .grid { display: grid; grid-template-columns: repeat(12, 1fr); gap: 16px; }
    .card { grid-column: span 6; background: #fff; border: 1px solid #d9e1e8; border-radius: 8px; padding: 16px; }
    .full { grid-column: 1 / -1; }
    table { width: 100%%; border-collapse: collapse; }
    th, td { text-align: left; padding: 9px 10px; border-bottom: 1px solid #d9e1e8; vertical-align: top; }
    th { color: #344054; background: #eef3f7; width: 260px; }
    img { width: 100%%; height: auto; border: 1px solid #d9e1e8; border-radius: 8px; background: white; }
    .metric-row { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin: 16px 0; }
    .metric { background: #fff; border: 1px solid #d9e1e8; border-radius: 8px; padding: 14px; }
    .metric strong { display: block; font-size: 24px; }
    .metric span { color: #667085; font-size: 13px; }
    @media (max-width: 900px) { .card { grid-column: 1 / -1; } .metric-row { grid-template-columns: 1fr 1fr; } }
  </style>
</head>
<body>
  <header>
    <h1>%s Seurat QC + PCA Report</h1>
    <p>%s · %s · naive top-variance HVGs after Seurat LogNormalize</p>
  </header>
  <main>
    <div class="metric-row">
      <div class="metric"><strong>%s</strong><span>cells before filtering</span></div>
      <div class="metric"><strong>%s</strong><span>features before filtering</span></div>
      <div class="metric"><strong>%s</strong><span>total UMIs/counts</span></div>
      <div class="metric"><strong>%s</strong><span>PC eigenvalue outliers</span></div>
    </div>
    <section class="grid">
      <article class="card full"><h2>Basic QC</h2>%s</article>
      <article class="card"><h2>UMI Rank Plot</h2><img src="%s" alt="UMI rank plot"></article>
      <article class="card"><h2>UMI Distribution</h2><img src="%s" alt="UMI distribution"></article>
      <article class="card"><h2>PCA Eigenvalues</h2><img src="%s" alt="PCA eigenvalue scree plot"></article>
      <article class="card"><h2>MP Bulk Fit</h2><img src="%s" alt="MP fit to PCA eigenvalue bulk"></article>
      <article class="card full"><h2>UMAP</h2><img src="%s" alt="UMAP"></article>
    </section>
  </main>
</body>
</html>',
  html_escape(dataset_name),
  html_escape(dataset_name),
  html_escape(technology),
  html_escape(assay_label),
  fmt_num(initial_cells),
  fmt_num(initial_features),
  fmt_num(sum(as.numeric(umi))),
  fmt_num(mp_fit$n_outliers_above_mp),
  table_html(qc_rows),
  rank_file,
  umi_file,
  scree_file,
  mp_file,
  umap_file
)

writeLines(html, file.path(out_dir, "index.html"))
message("Wrote report: ", file.path(out_dir, "index.html"))
