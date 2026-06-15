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

out_dir <- arg_value("--out-dir", "outputs/exploration/splatter_marker_depth_grid_even")
seed0 <- as.integer(arg_value("--seed", "20260608"))
p_eff <- as.integer(arg_value("--p-eff", "2000"))
n <- as.integer(arg_value("--n", "2000"))
rho <- as.numeric(arg_value("--rho", "0.5"))
de_fold <- as.numeric(arg_value("--de-fold", "4"))
max_plot_cells <- as.integer(arg_value("--max-plot-cells", "1600"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m_grid <- c(10L, 25L, 50L, 100L, 200L)
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

simulate_panel <- function(m, N, panel_id) {
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
  label <- ifelse(y == 1, "group 2", "group 1")
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
    m = m,
    N = N,
    cell = seq_len(n),
    group = label,
    kmeans = paste0("cluster ", km$cluster),
    PC1 = pc_scores[, 1],
    PC2 = pc_scores[, 2],
    seed = seed
  )

  if (nrow(plot_df) > max_plot_cells) {
    keep <- unlist(tapply(seq_len(nrow(plot_df)), plot_df$group, function(idx) {
      sample(idx, min(length(idx), ceiling(max_plot_cells / 2)))
    }))
    plot_df <- plot_df[sort(keep), ]
  }

  metrics <- data.frame(
    m = m,
    N = N,
    seed = seed,
    n = n,
    rho = rho,
    p_eff = p_eff,
    n_group2 = sum(y),
    mean_total_umi = mean(lib),
    scale_x = rho * m * N / p_eff,
    predicted_depth_scale = p_eff / (rho * m),
    max_auc = max(aucs, na.rm = TRUE),
    ari = ari,
    eig1 = svd_fit$d[1]^2 / max(1, n - 1),
    eig2 = svd_fit$d[2]^2 / max(1, n - 1),
    gap_1_2 = (svd_fit$d[1]^2) / (svd_fit$d[2]^2)
  )

  list(plot = plot_df, metrics = metrics)
}

grid <- expand.grid(m = m_grid, N = N_grid, KEEP.OUT.ATTRS = FALSE) %>%
  arrange(N, m) %>%
  mutate(panel_id = row_number())

message("Running ", nrow(grid), " marker-depth Splat panels.")
results <- lapply(seq_len(nrow(grid)), function(i) {
  simulate_panel(grid$m[i], grid$N[i], grid$panel_id[i])
})

plot_df <- bind_rows(lapply(results, `[[`, "plot"))
metrics <- bind_rows(lapply(results, `[[`, "metrics")) %>%
  mutate(
    m_label = paste0("m=", m),
    N_label = paste0("N=", N)
  )

plot_df <- plot_df %>%
  left_join(metrics %>% select(m, N, m_label, N_label), by = c("m", "N")) %>%
  mutate(
    group = factor(group, levels = c("group 1", "group 2")),
    kmeans = factor(kmeans),
    m_label = factor(m_label, levels = paste0("m=", m_grid)),
    N_label = factor(N_label, levels = paste0("N=", N_grid))
  )

ann_df <- metrics %>%
  mutate(
    m_label = factor(m_label, levels = paste0("m=", m_grid)),
    N_label = factor(N_label, levels = paste0("N=", N_grid)),
    label = paste0("AUC ", sprintf("%.2f", max_auc), "\nARI ", sprintf("%.2f", ari))
  )

write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
write.csv(plot_df, file.path(out_dir, "plot_points_sampled.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 11))

p_true <- ggplot(plot_df, aes(PC1, PC2)) +
  geom_point(
    data = subset(plot_df, group == "group 1"),
    color = "#9aa4b2",
    alpha = 0.35,
    size = 0.45
  ) +
  geom_point(
    data = subset(plot_df, group == "group 2"),
    color = "#d62728",
    alpha = 0.75,
    size = 0.55
  ) +
  geom_label(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 2.5,
    linewidth = 0,
    fill = "white",
    alpha = 0.75
  ) +
  facet_grid(N_label ~ m_label, scales = "free") +
  labs(
    title = "Balanced two-group Splat structure: depth by marker-count grid",
    subtitle = paste0(
      "Fixed p_eff=", p_eff,
      ", n=", n,
      ", rho=", rho,
      ", DE fold=", de_fold,
      "; columns vary marker genes, rows vary total UMI depth"
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
  geom_point(alpha = 0.55, size = 0.5) +
  geom_label(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 2.5,
    linewidth = 0,
    fill = "white",
    alpha = 0.75
  ) +
  facet_grid(N_label ~ m_label, scales = "free") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(
    title = "Same panels colored by k-means on PCs 1-10",
    subtitle = "ARI in each panel compares k-means clusters with true group labels",
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

p_heat <- ggplot(metrics, aes(factor(N), factor(m), fill = max_auc)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", max_auc)), size = 3) +
  scale_fill_viridis_c(limits = c(0.5, 1), option = "C") +
  labs(
    title = "Detectability over total UMI and marker-gene count",
    subtitle = "If the bound is inverse in m, the transition should move to lower N as m grows",
    x = "Expected total UMI per cell (N)",
    y = "Marker genes (m)",
    fill = "Max AUC"
  )

thresholds <- metrics %>%
  group_by(m, predicted_depth_scale) %>%
  arrange(N, .by_group = TRUE) %>%
  summarise(
    N_min_auc_0_8 = if (any(max_auc > 0.8)) min(N[max_auc > 0.8]) else NA_real_,
    N_min_ari_0_3 = if (any(ari > 0.3)) min(N[ari > 0.3]) else NA_real_,
    .groups = "drop"
  )
write.csv(thresholds, file.path(out_dir, "thresholds.csv"), row.names = FALSE)

p_threshold <- ggplot(thresholds, aes(predicted_depth_scale, N_min_auc_0_8)) +
  geom_point(size = 2.3, color = "#1f77b4") +
  geom_line(color = "#1f77b4") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Estimated N threshold versus p_eff / (rho m)",
    subtitle = "N_min is the first grid depth with max PC AUC > 0.8",
    x = "p_eff / (rho m)",
    y = "Estimated N_min"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(out_dir, "cluster_grid_true_labels.png"), p_true, width = 16, height = 12, dpi = 220)
ggsave(file.path(out_dir, "cluster_grid_kmeans.png"), p_kmeans, width = 16, height = 12, dpi = 220)
ggsave(file.path(out_dir, "auc_heatmap.png"), p_heat, width = 7.5, height = 4.8, dpi = 220)
ggsave(file.path(out_dir, "threshold_vs_marker_scale.png"), p_threshold, width = 6, height = 4.5, dpi = 220)

html_table <- function(x) {
  paste(capture.output(print(knitr::kable(
    x %>% mutate(across(where(is.numeric), ~ signif(.x, 4))),
    format = "html",
    table.attr = "class=\"data\""
  ))), collapse = "\n")
}

html <- paste0(
'<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Splatter marker-depth grid</title>
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
<h1>Splatter marker-depth grid</h1>
<p class="note">Balanced two-group structure; only total UMI depth and number of marker genes vary.</p>
<p>Fixed parameters: <code>p_eff=', p_eff, '</code>, <code>n=', n, '</code>, <code>rho=', rho, '</code>, <code>DE fold=', de_fold, '</code>. The background group has no differential expression; group 2 has marker probability <code>m / p_eff</code>.</p>

<h2>True group labels</h2>
<img src="cluster_grid_true_labels.png" alt="PCA grid colored by true labels">

<h2>K-means clusters</h2>
<img src="cluster_grid_kmeans.png" alt="PCA grid colored by k-means clusters">

<h2>Detection heatmap</h2>
<img src="auc_heatmap.png" alt="Max AUC heatmap">

<h2>Threshold scale</h2>
<img src="threshold_vs_marker_scale.png" alt="Threshold versus marker scale">
',
html_table(thresholds),
'<h2>Panel metrics</h2>',
html_table(metrics %>% select(m, N, mean_total_umi, scale_x, predicted_depth_scale, max_auc, ari, gap_1_2)),
'<p>Source files: <code>grid.csv</code>, <code>metrics.csv</code>, <code>thresholds.csv</code>, and <code>plot_points_sampled.csv</code>.</p>
</body>
</html>'
)

writeLines(html, file.path(out_dir, "splatter_marker_depth_grid.html"))
message("Wrote marker-depth grid outputs to ", out_dir)
