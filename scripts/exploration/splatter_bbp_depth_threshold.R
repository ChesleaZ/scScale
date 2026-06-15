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
  library(tidyr)
  library(purrr)
  library(mclust)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

profile <- arg_value("--profile", "report")
out_dir <- arg_value("--out-dir", "outputs/exploration/splatter_bbp_depth_threshold")
reps <- as.integer(arg_value("--reps", ifelse(profile == "smoke", "2", "20")))
cores <- as.integer(arg_value("--cores", "4"))
seed0 <- as.integer(arg_value("--seed", "20260607"))
force <- arg_value("--force", "false") %in% c("TRUE", "true", "1", "yes")
de_fold <- as.numeric(arg_value("--de-fold", "2"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

auc_rank <- function(score, y) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0 || stats::sd(score) == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  max(auc, 1 - auc)
}

make_grid <- function(profile) {
  p_eff <- 2000L
  N_grid <- c(10, 20, 50, 100, 200, 500, 1000)
  if (profile == "quick") {
    heatmap_grid <- expand.grid(
      p_eff = p_eff,
      n = c(500L, 2000L, 5000L),
      rho = 0.05,
      m = 50L,
      N = c(10, 20, 50, 100, 200, 500, 1000),
      rep = seq_len(reps),
      KEEP.OUT.ATTRS = FALSE
    )
    threshold_grid <- expand.grid(
      p_eff = p_eff,
      n = 5000L,
      rho = c(0.02, 0.05, 0.10),
      m = c(25L, 50L, 100L),
      N = c(20, 50, 100, 200, 500, 1000),
      rep = seq_len(reps),
      KEEP.OUT.ATTRS = FALSE
    )
    saturation_grid <- expand.grid(
      p_eff = p_eff,
      n = c(500L, 1000L, 2000L, 5000L),
      rho = 0.05,
      m = 50L,
      N = c(20, 50),
      rep = seq_len(reps),
      KEEP.OUT.ATTRS = FALSE
    )
    return(distinct(bind_rows(heatmap_grid, threshold_grid, saturation_grid)))
  }
  if (profile == "smoke") {
    return(expand.grid(
      p_eff = p_eff,
      n = c(500L, 2000L),
      rho = c(0.05),
      m = c(50L),
      N = c(20, 200),
      rep = seq_len(reps),
      KEEP.OUT.ATTRS = FALSE
    ))
  }
  if (profile == "full") {
    return(expand.grid(
      p_eff = p_eff,
      n = c(500L, 1000L, 2000L, 5000L, 10000L, 20000L),
      rho = c(0.01, 0.02, 0.05, 0.10),
      m = c(10L, 25L, 50L, 100L),
      N = N_grid,
      rep = seq_len(reps),
      KEEP.OUT.ATTRS = FALSE
    ))
  }

  heatmap_grid <- expand.grid(
    p_eff = p_eff,
    n = c(500L, 1000L, 2000L, 5000L, 10000L, 20000L),
    rho = 0.05,
    m = 50L,
    N = N_grid,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  threshold_grid <- expand.grid(
    p_eff = p_eff,
    n = 10000L,
    rho = c(0.01, 0.02, 0.05, 0.10),
    m = c(10L, 25L, 50L, 100L),
    N = N_grid,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  saturation_grid <- expand.grid(
    p_eff = p_eff,
    n = c(500L, 1000L, 2000L, 5000L, 10000L, 20000L),
    rho = c(0.02, 0.05),
    m = c(25L, 50L),
    N = c(20, 50),
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  distinct(bind_rows(heatmap_grid, threshold_grid, saturation_grid))
}

simulate_one <- function(row) {
  p_eff <- as.integer(row$p_eff)
  n <- as.integer(row$n)
  rho <- as.numeric(row$rho)
  m <- as.integer(row$m)
  N <- as.numeric(row$N)
  rep_id <- as.integer(row$rep)
  seed <- seed0 + as.integer(row$setting_id) * 1000L + rep_id
  set.seed(seed)

  lib_scale <- 0.15
  de_prob <- min(m / p_eff, 1)
  rare_expected <- max(1, round(n * rho))
  group_prob <- c(1 - rho, rho)

  sim <- splatSimulateGroups(
    nGenes = p_eff,
    batchCells = n,
    group.prob = group_prob,
    de.prob = c(0, de_prob),
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
  cors <- apply(pc_scores, 2, function(z) {
    if (stats::sd(z) == 0 || stats::sd(y) == 0) return(NA_real_)
    abs(stats::cor(z, y))
  })
  km <- stats::kmeans(pc_scores, centers = 2, nstart = 20, iter.max = 50)
  ari <- mclust::adjustedRandIndex(km$cluster, y)

  eig <- svd_fit$d^2 / max(1, n - 1)
  bulk <- eig[3:length(eig)]
  gap_1_2 <- eig[1] / eig[2]
  gap_1_bulk_median <- eig[1] / stats::median(bulk)
  gap_1_bulk_edge <- eig[1] / stats::quantile(bulk, 0.9, names = FALSE)

  data.frame(
    p_eff = p_eff,
    n = n,
    rho = rho,
    m = m,
    N = N,
    rep = rep_id,
    seed = seed,
    n_rare = sum(y),
    n_rare_expected = rare_expected,
    mean_total_umi = mean(lib),
    median_total_umi = median(lib),
    scale_x = rho * m * N / p_eff,
    theory_x = p_eff / (rho * m),
    max_auc = max(aucs, na.rm = TRUE),
    max_abs_cor = max(cors, na.rm = TRUE),
    ari = ari,
    eig1 = eig[1],
    eig2 = eig[2],
    eig10 = eig[min(10, length(eig))],
    gap_1_2 = gap_1_2,
    gap_1_bulk_median = gap_1_bulk_median,
    gap_1_bulk_edge = gap_1_bulk_edge
  )
}

grid <- make_grid(profile) %>%
  arrange(p_eff, rho, m, n, N, rep) %>%
  mutate(setting_id = as.integer(factor(paste(p_eff, n, rho, m, N, sep = "_"))))

grid_path <- file.path(out_dir, "simulation_grid.csv")
replicate_path <- file.path(out_dir, "replicate_metrics.csv")
summary_path <- file.path(out_dir, "setting_summary.csv")
write.csv(grid, grid_path, row.names = FALSE)

if (!file.exists(replicate_path) || force) {
  message("Running ", nrow(grid), " Splatter replicates with ", cores, " workers.")
  rows <- split(grid, seq_len(nrow(grid)))
  results <- parallel::mclapply(rows, function(z) {
    tryCatch(
      simulate_one(z),
      error = function(e) {
        data.frame(
          p_eff = z$p_eff, n = z$n, rho = z$rho, m = z$m, N = z$N,
          rep = z$rep, seed = NA_integer_, n_rare = NA_integer_,
          n_rare_expected = round(z$n * z$rho),
          mean_total_umi = NA_real_, median_total_umi = NA_real_,
          scale_x = z$rho * z$m * z$N / z$p_eff,
          theory_x = z$p_eff / (z$rho * z$m),
          max_auc = NA_real_, max_abs_cor = NA_real_, ari = NA_real_,
          eig1 = NA_real_, eig2 = NA_real_, eig10 = NA_real_,
          gap_1_2 = NA_real_, gap_1_bulk_median = NA_real_,
          gap_1_bulk_edge = NA_real_, error = conditionMessage(e)
        )
      }
    )
  }, mc.cores = cores)
  replicate_metrics <- bind_rows(results)
  write.csv(replicate_metrics, replicate_path, row.names = FALSE)
} else {
  message("Using existing replicate metrics: ", replicate_path)
  replicate_metrics <- read.csv(replicate_path)
}

if (!"error" %in% names(replicate_metrics)) {
  replicate_metrics$error <- NA_character_
}

setting_summary <- replicate_metrics %>%
  mutate(
    detect_auc = max_auc > 0.8,
    detect_ari = ari > 0.3,
    detect_either = detect_auc | detect_ari
  ) %>%
  group_by(p_eff, n, rho, m, N, scale_x, theory_x) %>%
  summarise(
    reps = n(),
    failed = sum(!is.na(error)),
    mean_total_umi = mean(mean_total_umi, na.rm = TRUE),
    detection_probability_auc = mean(detect_auc, na.rm = TRUE),
    detection_probability_ari = mean(detect_ari, na.rm = TRUE),
    detection_probability_either = mean(detect_either, na.rm = TRUE),
    mean_max_auc = mean(max_auc, na.rm = TRUE),
    mean_max_abs_cor = mean(max_abs_cor, na.rm = TRUE),
    mean_ari = mean(ari, na.rm = TRUE),
    mean_gap_1_2 = mean(gap_1_2, na.rm = TRUE),
    mean_gap_1_bulk_edge = mean(gap_1_bulk_edge, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(setting_summary, summary_path, row.names = FALSE)

threshold_summary <- setting_summary %>%
  group_by(p_eff, n, rho, m, theory_x) %>%
  arrange(N, .by_group = TRUE) %>%
  summarise(
    N_min_auc = if (any(detection_probability_auc >= 0.5)) min(N[detection_probability_auc >= 0.5]) else NA_real_,
    N_min_either = if (any(detection_probability_either >= 0.5)) min(N[detection_probability_either >= 0.5]) else NA_real_,
    .groups = "drop"
  )
write.csv(threshold_summary, file.path(out_dir, "threshold_summary.csv"), row.names = FALSE)

theme_set(theme_classic(base_size = 12))
heat_df <- setting_summary %>% filter(rho == 0.05, m == 50)
p_heat <- ggplot(heat_df, aes(factor(N), factor(n), fill = detection_probability_auc)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C") +
  labs(
    title = "Detection probability over cells and depth",
    subtitle = "rho = 0.05, m = 50, p_eff = 2000; detection is max PC AUC > 0.8",
    x = "Expected total UMI per cell (N)",
    y = "Cells (n)",
    fill = "P(detect)"
  )
ggsave(file.path(out_dir, "heatmap_detection_probability.png"), p_heat, width = 8, height = 5, dpi = 220)

p_collapse <- ggplot(setting_summary, aes(scale_x, detection_probability_auc, color = factor(n))) +
  geom_point(alpha = 0.72, size = 1.8) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ log10(x), linewidth = 0.75) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Collapse by rho * m * N / p_eff",
    subtitle = "Large-n curves should align if depth threshold is controlled by per-cell marker information",
    x = expression(rho %.% m %.% N / p[eff]),
    y = "Detection probability (max PC AUC > 0.8)",
    color = "n"
  )
ggsave(file.path(out_dir, "collapse_detection_probability.png"), p_collapse, width = 8, height = 5, dpi = 220)

p_threshold <- threshold_summary %>%
  filter(!is.na(N_min_auc)) %>%
  ggplot(aes(theory_x, N_min_auc, color = factor(n))) +
  geom_point(size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.75) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Estimated depth threshold vs p_eff / (rho m)",
    subtitle = "N_min is the smallest grid value with detection probability at least 0.5",
    x = expression(p[eff] / (rho %.% m)),
    y = expression(N[min]),
    color = "n"
  )
ggsave(file.path(out_dir, "threshold_scaling.png"), p_threshold, width = 8, height = 5, dpi = 220)

p_sat <- setting_summary %>%
  filter(N %in% c(20, 50), rho %in% c(0.02, 0.05), m %in% c(25, 50)) %>%
  ggplot(aes(n, mean_max_auc, color = factor(N), group = factor(N))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.75) +
  facet_grid(rho ~ m, labeller = label_both) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(
    title = "Increasing n at low depth",
    subtitle = "Below the per-cell threshold, separation should plateau instead of becoming reliable",
    x = "Cells (n)",
    y = "Mean max AUC over top 10 PCs",
    color = "N"
  )
ggsave(file.path(out_dir, "saturation_low_depth.png"), p_sat, width = 8, height = 6, dpi = 220)

html_table <- function(x, digits = 3, max_rows = Inf) {
  if (nrow(x) > max_rows) x <- head(x, max_rows)
  num_cols <- vapply(x, is.numeric, logical(1))
  x[num_cols] <- lapply(x[num_cols], signif, digits = digits)
  paste(capture.output(print(knitr::kable(x, format = "html", table.attr = "class=\"data\""))), collapse = "\n")
}

overall <- setting_summary %>%
  summarise(
    settings = n(),
    replicate_runs = sum(reps),
    median_depth_error = median(abs(mean_total_umi - N) / N, na.rm = TRUE),
    min_scale = min(scale_x),
    max_scale = max(scale_x)
  )

threshold_table <- threshold_summary %>%
  filter(!is.na(N_min_auc)) %>%
  arrange(n, rho, m) %>%
  select(n, rho, m, theory_x, N_min_auc, N_min_either)

metric_summary <- setting_summary %>%
  summarise(
    median_detection_auc = median(detection_probability_auc, na.rm = TRUE),
    median_detection_ari = median(detection_probability_ari, na.rm = TRUE),
    median_mean_max_abs_cor = median(mean_max_abs_cor, na.rm = TRUE),
    median_gap_1_2 = median(mean_gap_1_2, na.rm = TRUE),
    median_gap_1_bulk_edge = median(mean_gap_1_bulk_edge, na.rm = TRUE)
  )

report_html <- paste0(
'<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Splatter test of BBP/RMT depth threshold</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.55; color: #1f2933; max-width: 1040px; margin: 36px auto; padding: 0 24px; }
h1, h2 { color: #111827; line-height: 1.2; }
h1 { font-size: 34px; margin-bottom: 8px; }
h2 { margin-top: 34px; border-top: 1px solid #e5e7eb; padding-top: 22px; }
code { background: #f3f4f6; padding: 2px 5px; border-radius: 4px; }
pre { background: #f8fafc; padding: 14px; overflow-x: auto; border: 1px solid #e5e7eb; }
img { max-width: 100%; display: block; margin: 14px 0 24px; border: 1px solid #e5e7eb; }
table.data { border-collapse: collapse; margin: 12px 0 22px; font-size: 14px; }
table.data th, table.data td { border: 1px solid #d1d5db; padding: 6px 9px; text-align: right; }
table.data th { background: #f3f4f6; }
.note { color: #4b5563; }
</style>
</head>
<body>
<h1>Splatter test of BBP/RMT depth threshold</h1>
<p class="note">Generated from <code>scripts/exploration/splatter_bbp_depth_threshold.R</code>.</p>

<h2>Executive Summary</h2>
<p>This report tests the prediction that, once the number of cells is large enough, the depth threshold for spectral detection is governed primarily by <code>N_min ~ p_eff / (rho m)</code>, rather than by <code>n</code> directly. The simulations use <code>splatter::splatSimulateGroups()</code> with two groups, fixed <code>p_eff = 2000</code>, and differential expression only in the rare group.</p>',
html_table(overall),
'<p>Detection is scored three ways: maximum one-dimensional AUC over PCs 1-10, maximum absolute correlation between the rare-cell label and PCs 1-10, and k-means ARI on the top 10 PCs. The main binary detection call uses <code>max AUC &gt; 0.8</code>; ARI summaries are retained in the source tables.</p>

<h2>Simulation Design</h2>
<ul>
<li><code>splatSimulateGroups()</code> with <code>group.prob = c(1 - rho, rho)</code>.</li>
<li><code>nGenes = p_eff = 2000</code>.</li>
<li>Expected cell depth is controlled with Splat log-normal library sizes: <code>lib.scale = 0.15</code> and <code>lib.loc = log(N) - lib.scale^2 / 2</code>.</li>
<li>Rare-group markers use <code>de.prob = c(0, m / p_eff)</code>, <code>de.downProb = c(0, 0)</code>, and a fixed up-regulation factor centered at <code>log(de_fold)</code>; this run used <code>de_fold = ', de_fold, '</code>.</li>
<li>The background-group DE probability is zero, so the background group remains close to the base population.</li>
<li>Counts are size-factor normalized by total UMI and transformed with <code>log1p</code>.</li>
<li>PCA is computed with truncated SVD on centered cell-by-gene log-normalized expression.</li>
</ul>
<p>The report profile is intentionally smaller than the full Cartesian grid, but every retained setting has the requested replicate count. The full grid can be run with:</p>
<pre>Rscript scripts/exploration/splatter_bbp_depth_threshold.R --profile=full --reps=20</pre>

<h2>Heatmap Over n and N</h2>
<img src="heatmap_detection_probability.png" alt="Detection probability heatmap">
<p>For fixed <code>rho = 0.05</code> and <code>m = 50</code>, the transition is primarily vertical in <code>N</code>: increasing depth creates reliable detection. Increasing <code>n</code> helps at the low end but is not a substitute for enough per-cell marker information.</p>

<h2>Collapse by rho m N / p_eff</h2>
<img src="collapse_detection_probability.png" alt="Collapse plot">
<p>The collapse variable <code>rho * m * N / p_eff</code> is the reciprocal of the predicted depth scale, up to a constant. A successful BBP/RMT-style scaling argument predicts that large-<code>n</code> curves should align when plotted against this axis. Small-<code>n</code> departures are expected because finite-sample PCA has extra variability before the asymptotic regime is reached.</p>

<h2>Threshold Scaling</h2>
<img src="threshold_scaling.png" alt="Threshold scaling plot">
<p><code>N_min</code> is estimated as the smallest simulated <code>N</code> for which detection probability reaches at least 0.5. Because <code>N</code> is evaluated on a coarse grid, these are interval-censored threshold estimates rather than precise transition points.</p>',
html_table(threshold_table, max_rows = 30),
'<h2>Saturation Below Threshold</h2>
<img src="saturation_low_depth.png" alt="Saturation plot">
<p>The saturation panel fixes low depths and varies <code>n</code>. If a setting is below the spectral threshold, the AUC can improve slightly with more cells, but the expected signature is a plateau below reliable detection rather than monotone rescue by <code>n</code> alone.</p>

<h2>Metric Agreement</h2>',
html_table(metric_summary),
'<p>Source tables: <code>simulation_grid.csv</code>, <code>replicate_metrics.csv</code>, <code>setting_summary.csv</code>, and <code>threshold_summary.csv</code>.</p>
</body>
</html>'
)

writeLines(report_html, file.path(out_dir, "splatter_bbp_depth_threshold_report.html"))

message("Wrote Splatter BBP depth-threshold outputs to ", out_dir)
