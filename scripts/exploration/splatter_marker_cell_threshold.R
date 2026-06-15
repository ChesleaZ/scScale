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

out_dir <- arg_value("--out-dir", "outputs/exploration/splatter_marker_cell_threshold_even")
seed0 <- as.integer(arg_value("--seed", "20260608"))
p_eff <- as.integer(arg_value("--p-eff", "2000"))
N <- as.numeric(arg_value("--N", "500"))
rho <- as.numeric(arg_value("--rho", "0.5"))
de_fold <- as.numeric(arg_value("--de-fold", "4"))
reps <- as.integer(arg_value("--reps", "3"))
max_plot_cells <- as.integer(arg_value("--max-plot-cells", "1600"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m_grid <- c(10L, 25L, 50L, 100L, 200L)
n_grid <- c(100L, 200L, 500L, 1000L, 2000L, 5000L)

auc_rank <- function(score, y) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0 || stats::sd(score) == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  max(auc, 1 - auc)
}

simulate_one <- function(n, m, rep_id, panel_id) {
  seed <- seed0 + panel_id * 100L + rep_id
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
  km <- stats::kmeans(pc_scores, centers = 2, nstart = 25, iter.max = 80)
  ari <- mclust::adjustedRandIndex(km$cluster, y)

  plot_df <- NULL
  if (rep_id == 1) {
    plot_df <- data.frame(
      n = n,
      m = m,
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
  }

  metrics <- data.frame(
    n = n,
    m = m,
    rep = rep_id,
    seed = seed,
    rho = rho,
    N = N,
    p_eff = p_eff,
    n_group2 = sum(y),
    mean_total_umi = mean(lib),
    marker_scale = rho * m * N / p_eff,
    inverse_marker_scale = p_eff / (rho * m),
    max_auc = max(aucs, na.rm = TRUE),
    ari = ari,
    eig1 = svd_fit$d[1]^2 / max(1, n - 1),
    eig2 = svd_fit$d[2]^2 / max(1, n - 1),
    gap_1_2 = (svd_fit$d[1]^2) / (svd_fit$d[2]^2)
  )

  list(plot = plot_df, metrics = metrics)
}

grid <- expand.grid(n = n_grid, m = m_grid, rep = seq_len(reps), KEEP.OUT.ATTRS = FALSE) %>%
  arrange(m, n, rep) %>%
  mutate(panel_id = as.integer(factor(paste(n, m, sep = "_"))))

message("Running ", nrow(grid), " marker-cell Splat replicates.")
results <- lapply(seq_len(nrow(grid)), function(i) {
  simulate_one(grid$n[i], grid$m[i], grid$rep[i], grid$panel_id[i])
})

metrics <- bind_rows(lapply(results, `[[`, "metrics"))
plot_df <- bind_rows(lapply(results, `[[`, "plot"))

summary_df <- metrics %>%
  group_by(m, n, inverse_marker_scale) %>%
  summarise(
    reps = n(),
    mean_ari = mean(ari),
    median_ari = median(ari),
    p_ari_0_9 = mean(ari >= 0.9),
    mean_max_auc = mean(max_auc),
    median_max_auc = median(max_auc),
    .groups = "drop"
  )

thresholds <- summary_df %>%
  group_by(m, inverse_marker_scale) %>%
  arrange(n, .by_group = TRUE) %>%
  summarise(
    n_min_median_ari_0_9 = if (any(median_ari >= 0.9)) min(n[median_ari >= 0.9]) else NA_real_,
    n_min_p_ari_0_9_ge_half = if (any(p_ari_0_9 >= 0.5)) min(n[p_ari_0_9 >= 0.5]) else NA_real_,
    .groups = "drop"
  )

metrics_for_labels <- metrics %>%
  filter(rep == 1) %>%
  mutate(
    m_label = paste0("m=", m),
    n_label = paste0("n=", n),
    label = paste0("ARI ", sprintf("%.2f", ari), "\nAUC ", sprintf("%.2f", max_auc))
  )

plot_df <- plot_df %>%
  mutate(
    group = factor(group, levels = c("group 1", "group 2")),
    kmeans = factor(kmeans),
    m_label = factor(paste0("m=", m), levels = paste0("m=", m_grid)),
    n_label = factor(paste0("n=", n), levels = paste0("n=", n_grid))
  )

ann_df <- metrics_for_labels %>%
  mutate(
    m_label = factor(m_label, levels = paste0("m=", m_grid)),
    n_label = factor(n_label, levels = paste0("n=", n_grid))
  )

write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
write.csv(metrics, file.path(out_dir, "replicate_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "setting_summary.csv"), row.names = FALSE)
write.csv(thresholds, file.path(out_dir, "thresholds.csv"), row.names = FALSE)
write.csv(plot_df, file.path(out_dir, "plot_points_sampled_rep1.csv"), row.names = FALSE)

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
  facet_grid(n_label ~ m_label, scales = "free") +
  labs(
    title = "Balanced two-group Splat structure: cell count by marker-count grid",
    subtitle = paste0(
      "Fixed p_eff=", p_eff,
      ", N=", N,
      ", rho=", rho,
      ", DE fold=", de_fold,
      "; columns vary marker genes, rows vary cell count"
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
  facet_grid(n_label ~ m_label, scales = "free") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(
    title = "Same panels colored by k-means on PCs 1-10",
    subtitle = "ARI in each panel is for the displayed replicate; thresholds use all replicates",
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

p_heat <- ggplot(summary_df, aes(factor(n), factor(m), fill = median_ari)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", median_ari)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C") +
  labs(
    title = "Median ARI over replicates",
    subtitle = "Detection threshold: median ARI >= 0.9",
    x = "Cells (n)",
    y = "Marker genes (m)",
    fill = "Median ARI"
  )

p_threshold <- ggplot(thresholds, aes(m, n_min_median_ari_0_9)) +
  geom_point(size = 2.3, color = "#1f77b4") +
  geom_line(color = "#1f77b4") +
  scale_x_log10(breaks = m_grid) +
  scale_y_log10(breaks = n_grid) +
  labs(
    title = "Minimum n for median ARI >= 0.9",
    subtitle = paste0("Fixed N=", N, "; inverse-marker trend should slope downward if more markers reduce sample need"),
    x = "Marker genes (m)",
    y = "Minimum n"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(out_dir, "cluster_grid_true_labels.png"), p_true, width = 16, height = 12, dpi = 220)
ggsave(file.path(out_dir, "cluster_grid_kmeans.png"), p_kmeans, width = 16, height = 12, dpi = 220)
ggsave(file.path(out_dir, "median_ari_heatmap.png"), p_heat, width = 7.5, height = 4.8, dpi = 220)
ggsave(file.path(out_dir, "n_threshold_vs_marker_count.png"), p_threshold, width = 6, height = 4.5, dpi = 220)

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
<title>Splatter marker-cell ARI threshold</title>
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
<h1>Minimum n for ARI >= 0.9 by marker count</h1>
<p class="note">Balanced two-group structure; fixed total UMI depth; cell count and number of marker genes vary.</p>
<p>Fixed parameters: <code>p_eff=', p_eff, '</code>, <code>N=', N, '</code>, <code>rho=', rho, '</code>, <code>DE fold=', de_fold, '</code>, <code>reps=', reps, '</code>. The background group has no differential expression; group 2 has marker probability <code>m / p_eff</code>.</p>

<h2>True group labels, replicate 1</h2>
<img src="cluster_grid_true_labels.png" alt="PCA grid colored by true labels">

<h2>K-means clusters, replicate 1</h2>
<img src="cluster_grid_kmeans.png" alt="PCA grid colored by k-means clusters">

<h2>Median ARI heatmap</h2>
<img src="median_ari_heatmap.png" alt="Median ARI heatmap">

<h2>Threshold trend</h2>
<img src="n_threshold_vs_marker_count.png" alt="Minimum n threshold versus marker count">
',
html_table(thresholds),
'<h2>Setting summary</h2>',
html_table(summary_df),
'<p>Source files: <code>replicate_metrics.csv</code>, <code>setting_summary.csv</code>, <code>thresholds.csv</code>, and <code>plot_points_sampled_rep1.csv</code>.</p>
</body>
</html>'
)

writeLines(html, file.path(out_dir, "splatter_marker_cell_threshold.html"))
message("Wrote marker-cell threshold outputs to ", out_dir)
