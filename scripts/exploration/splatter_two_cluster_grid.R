local_lib <- normalizePath(".rlib/R", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(splatter)
  library(SummarizedExperiment)
  library(irlba)
  library(ggplot2)
  library(dplyr)
  library(mclust)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

out_dir <- arg_value("--out-dir", "outputs/exploration/splatter_two_cluster_grid")
seed0 <- as.integer(arg_value("--seed", "20260608"))
p_eff <- as.integer(arg_value("--p-eff", "2000"))
rho <- as.numeric(arg_value("--rho", "0.05"))
m <- as.integer(arg_value("--m", "50"))
de_fold <- as.numeric(arg_value("--de-fold", "4"))
max_plot_cells <- as.integer(arg_value("--max-plot-cells", "1200"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_grid <- c(500L, 1000L, 2000L, 5000L)
N_grid <- c(20, 50, 100, 200, 500, 1000)

auc_rank <- function(score, y) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0 || stats::sd(score) == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  max(auc, 1 - auc)
}

simulate_panel <- function(n, N, panel_id) {
  seed <- seed0 + panel_id
  set.seed(seed)

  lib_scale <- 0.15
  sim <- splatSimulateGroups(
    nGenes = p_eff,
    batchCells = n,
    group.prob = c(1 - rho, rho),
    de.prob = c(0, min(m / p_eff, 1)),
    de.downProb = c(0, 0),
    de.facLoc = c(0, log(de_fold)),
    de.facScale = c(0, 0.1),
    lib.loc = log(N) - lib_scale^2 / 2,
    lib.scale = lib_scale,
    out.prob = 0,
    dropout.type = "none",
    bcv.common = 0.1,
    seed = seed,
    verbose = FALSE
  )

  counts <- assay(sim, "counts")
  y <- as.integer(colData(sim)$Group == "Group2")
  label <- ifelse(y == 1, "rare", "background")
  lib <- Matrix::colSums(counts)
  median_lib <- stats::median(lib[lib > 0])
  sf <- lib / median_lib
  sf[!is.finite(sf) | sf <= 0] <- 1
  x <- log1p(t(t(as.matrix(counts)) / sf))
  x_cell_gene <- scale(t(x), center = TRUE, scale = FALSE)

  svd_fit <- irlba(x_cell_gene, nv = 10, nu = 10)
  pc_scores <- svd_fit$u %*% diag(svd_fit$d, nrow = length(svd_fit$d))
  colnames(pc_scores) <- paste0("PC", seq_len(ncol(pc_scores)))

  aucs <- apply(pc_scores, 2, auc_rank, y = y)
  km <- stats::kmeans(pc_scores, centers = 2, nstart = 25, iter.max = 60)
  ari <- mclust::adjustedRandIndex(km$cluster, y)

  plot_df <- data.frame(
    n = n,
    N = N,
    cell = seq_len(n),
    group = label,
    kmeans = paste0("cluster ", km$cluster),
    PC1 = pc_scores[, 1],
    PC2 = pc_scores[, 2],
    seed = seed
  )

  rare_idx <- which(plot_df$group == "rare")
  bg_idx <- which(plot_df$group == "background")
  bg_keep <- sample(bg_idx, min(length(bg_idx), max_plot_cells - min(length(rare_idx), 300)))
  rare_keep <- if (length(rare_idx) > 300) sample(rare_idx, 300) else rare_idx
  keep <- sort(unique(c(bg_keep, rare_keep)))
  plot_df <- plot_df[keep, ]

  metrics <- data.frame(
    n = n,
    N = N,
    seed = seed,
    n_rare = sum(y),
    mean_total_umi = mean(lib),
    scale_x = rho * m * N / p_eff,
    max_auc = max(aucs, na.rm = TRUE),
    ari = ari,
    eig1 = svd_fit$d[1]^2 / max(1, n - 1),
    eig2 = svd_fit$d[2]^2 / max(1, n - 1),
    gap_1_2 = (svd_fit$d[1]^2) / (svd_fit$d[2]^2)
  )

  list(plot = plot_df, metrics = metrics)
}

grid <- expand.grid(n = n_grid, N = N_grid, KEEP.OUT.ATTRS = FALSE) %>%
  arrange(N, n) %>%
  mutate(panel_id = row_number())

message("Running ", nrow(grid), " fixed-structure Splat panels.")
results <- lapply(seq_len(nrow(grid)), function(i) {
  simulate_panel(grid$n[i], grid$N[i], grid$panel_id[i])
})

plot_df <- bind_rows(lapply(results, `[[`, "plot"))
metrics <- bind_rows(lapply(results, `[[`, "metrics")) %>%
  mutate(
    n_label = paste0("n=", n),
    N_label = paste0("N=", N)
  )

plot_df <- plot_df %>%
  left_join(metrics %>% select(n, N, n_label, N_label), by = c("n", "N")) %>%
  mutate(
    group = factor(group, levels = c("background", "rare")),
    kmeans = factor(kmeans),
    n_label = factor(n_label, levels = paste0("n=", n_grid)),
    N_label = factor(N_label, levels = paste0("N=", N_grid))
  )

ann_df <- metrics %>%
  mutate(
    n_label = factor(n_label, levels = paste0("n=", n_grid)),
    N_label = factor(N_label, levels = paste0("N=", N_grid)),
    label = paste0("AUC ", sprintf("%.2f", max_auc), "\nARI ", sprintf("%.2f", ari))
  )

write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
write.csv(plot_df, file.path(out_dir, "plot_points_sampled.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 11))

p_true <- ggplot(plot_df, aes(PC1, PC2)) +
  geom_point(
    data = subset(plot_df, group == "background"),
    color = "#9aa4b2",
    alpha = 0.35,
    size = 0.45
  ) +
  geom_point(
    data = subset(plot_df, group == "rare"),
    color = "#d62728",
    alpha = 0.9,
    size = 0.8
  ) +
  geom_label(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 2.5,
    label.size = 0,
    fill = "white",
    alpha = 0.75
  ) +
  facet_grid(N_label ~ n_label, scales = "free") +
  labs(
    title = "Two-group Splat structure: PCA cluster grid",
    subtitle = paste0(
      "Fixed p_eff=", p_eff,
      ", rho=", rho,
      ", m=", m,
      ", DE fold=", de_fold,
      "; columns vary cell count, rows vary total UMI depth"
    ),
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_kmeans <- ggplot(plot_df, aes(PC1, PC2, color = kmeans)) +
  geom_point(alpha = 0.55, size = 0.55) +
  geom_label(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 2.5,
    label.size = 0,
    fill = "white",
    alpha = 0.75
  ) +
  facet_grid(N_label ~ n_label, scales = "free") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(
    title = "Same panels colored by k-means on PCs 1-10",
    subtitle = "ARI in each column header compares k-means clusters with the true rare/background labels",
    x = "PC1",
    y = "PC2",
    color = NULL
  ) +
  theme(
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  )

p_heat <- ggplot(metrics, aes(factor(N), factor(n), fill = max_auc)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", max_auc)), size = 3) +
  scale_fill_viridis_c(limits = c(0.5, 1), option = "C") +
  labs(
    title = "PC detectability for the same fixed two-group structure",
    x = "Expected total UMI per cell (N)",
    y = "Cells (n)",
    fill = "Max AUC"
  )

ggsave(file.path(out_dir, "cluster_grid_true_labels.png"), p_true, width = 15, height = 12, dpi = 220)
ggsave(file.path(out_dir, "cluster_grid_kmeans.png"), p_kmeans, width = 15, height = 12, dpi = 220)
ggsave(file.path(out_dir, "auc_heatmap.png"), p_heat, width = 7, height = 4.5, dpi = 220)

metrics_table <- paste(
  capture.output(print(knitr::kable(
    metrics %>%
      arrange(n, N) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4))),
    format = "html",
    table.attr = "class=\"data\""
  ))),
  collapse = "\n"
)

html <- paste0(
'<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Splatter two-cluster grid</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.55; color: #1f2933; max-width: 1280px; margin: 32px auto; padding: 0 24px; }
h1, h2 { color: #111827; line-height: 1.2; }
h2 { border-top: 1px solid #e5e7eb; padding-top: 22px; margin-top: 34px; }
code { background: #f3f4f6; padding: 2px 5px; border-radius: 4px; }
img { max-width: 100%; display: block; margin: 12px 0 24px; border: 1px solid #e5e7eb; }
table.data { border-collapse: collapse; margin: 12px 0 22px; font-size: 13px; }
table.data th, table.data td { border: 1px solid #d1d5db; padding: 5px 8px; text-align: right; }
table.data th { background: #f3f4f6; }
.note { color: #4b5563; }
</style>
</head>
<body>
<h1>Splatter two-cluster grid over N and n</h1>
<p class="note">One fixed two-group structure; only total UMI depth and cell count vary.</p>
<p>Fixed parameters: <code>p_eff=', p_eff, '</code>, <code>rho=', rho, '</code>, <code>m=', m, '</code>, <code>DE fold=', de_fold, '</code>. The background group has no differential expression; the rare group has marker probability <code>m / p_eff</code>. Each panel is one Splat simulation seed.</p>

<h2>True rare/background labels</h2>
<img src="cluster_grid_true_labels.png" alt="PCA grid colored by true labels">

<h2>K-means clusters</h2>
<img src="cluster_grid_kmeans.png" alt="PCA grid colored by k-means clusters">

<h2>Detection heatmap</h2>
<img src="auc_heatmap.png" alt="Max AUC heatmap">

<h2>Metrics</h2>',
metrics_table,
'<p>Source files: <code>grid.csv</code>, <code>metrics.csv</code>, and <code>plot_points_sampled.csv</code>.</p>
</body>
</html>'
)

writeLines(html, file.path(out_dir, "splatter_two_cluster_grid.html"))
message("Wrote two-cluster grid outputs to ", out_dir)
