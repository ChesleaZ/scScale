library(ggplot2)

out_dir <- "outputs/exploration/spike_scaling_by_dataset"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- data.frame(
  dataset = c(
    "Jurkat",
    "K562",
    "Norman-Weissman",
    "Replogle RPE1"
  ),
  eigen_file = c(
    "outputs/exploration/jurkat_plain_hvg_by_cells/jurkat_plain_hvg_by_cells_plain_hvg_eigenvalues_by_cells.csv",
    "outputs/exploration/k562_plain_hvg_by_cells/k562_plain_hvg_by_cells_plain_hvg_eigenvalues_by_cells.csv",
    "outputs/exploration/norman_weissman_2019_plain_hvg_by_cells/norman_weissman_2019_plain_hvg_by_cells_hvg_eigenvalues_by_cells.csv",
    "outputs/exploration/replogle_rpe1_plain_hvg_by_cells/replogle_rpe1_plain_hvg_by_cells_hvg_eigenvalues_by_cells.csv"
  ),
  mp_file = c(
    "outputs/exploration/jurkat_plain_hvg_mp_by_ncells/jurkat_mp_by_ncells_parameters.csv",
    "outputs/exploration/k562_plain_hvg_mp_by_ncells/jurkat_mp_by_ncells_parameters.csv",
    "outputs/exploration/norman_weissman_2019_plain_hvg_mp_by_ncells/jurkat_mp_by_ncells_parameters.csv",
    "outputs/exploration/replogle_rpe1_plain_hvg_mp_by_ncells/jurkat_mp_by_ncells_parameters.csv"
  ),
  stringsAsFactors = FALSE
)

count_spikes <- function(dataset, eigen_file, mp_file) {
  eig <- read.csv(eigen_file)
  fit <- read.csv(mp_file)

  rows <- vector("list", nrow(fit))
  for (i in seq_len(nrow(fit))) {
    fit_row <- fit[i, ]
    vals <- eig$eigenvalue[eig$n_cells == fit_row$n_cells]
    vals <- vals[is.finite(vals)]
    positive_vals <- vals[vals > 1e-10]
    spike_count <- sum(vals > fit_row$mp_upper)

    rows[[i]] <- data.frame(
      dataset = dataset,
      n_cells = fit_row$n_cells,
      p = fit_row$p,
      fit_type = fit_row$fit_type,
      gamma = fit_row$gamma,
      gamma_theory = fit_row$gamma_theory,
      sigma2 = fit_row$sigma2,
      mp_upper = fit_row$mp_upper,
      n_eigenvalues = length(vals),
      n_positive = length(positive_vals),
      spike_count = spike_count,
      spike_fraction_all = spike_count / length(vals),
      spike_fraction_positive = spike_count / length(positive_vals),
      max_eigenvalue = max(vals),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

spike_df <- do.call(
  rbind,
  lapply(seq_len(nrow(datasets)), function(i) {
    count_spikes(datasets$dataset[i], datasets$eigen_file[i], datasets$mp_file[i])
  })
)

spike_df$dataset <- factor(spike_df$dataset, levels = datasets$dataset)
spike_df$fit_type <- factor(
  spike_df$fit_type,
  levels = c("fixed_gamma", "free_gamma"),
  labels = c("Fixed gamma = p / n", "Free gamma")
)
spike_df$spike_count_for_plot <- spike_df$spike_count
spike_df$spike_fraction_all_for_plot <- spike_df$spike_fraction_all
fixed_gamma_singular <- spike_df$fit_type == "Fixed gamma = p / n" &
  abs(spike_df$gamma_theory - 1) < 1e-8
spike_df$spike_count_for_plot[fixed_gamma_singular] <- NA
spike_df$spike_fraction_all_for_plot[fixed_gamma_singular] <- NA

write.csv(
  spike_df,
  file.path(out_dir, "spike_count_summary.csv"),
  row.names = FALSE
)

theme_set(theme_classic(base_size = 13))

count_plot <- ggplot(spike_df, aes(n_cells, spike_count_for_plot, color = dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ fit_type, nrow = 1) +
  labs(
    title = "Eigenvalue spikes above the MP upper edge",
    subtitle = "Spike = eigenvalue greater than the fitted MP upper support",
    x = "Number of cells",
    y = "Spike count",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
)

ggsave(
  file.path(out_dir, "spike_count_vs_ncells.png"),
  count_plot,
  width = 10,
  height = 5.8,
  dpi = 240
)

high_n_plot <- ggplot(
  spike_df[spike_df$n_cells >= 5000, ],
  aes(n_cells, spike_count_for_plot, color = dataset)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ fit_type, nrow = 1) +
  labs(
    title = "Eigenvalue spikes above the MP upper edge, n >= 5,000",
    subtitle = "Low-n rank edge cases removed from view",
    x = "Number of cells",
    y = "Spike count",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(out_dir, "spike_count_vs_ncells_high_n.png"),
  high_n_plot,
  width = 10,
  height = 5.8,
  dpi = 240
)

fraction_plot <- ggplot(spike_df, aes(n_cells, spike_fraction_all_for_plot, color = dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  facet_wrap(~ fit_type, nrow = 1) +
  labs(
    title = "Fraction of eigenvalues above the MP upper edge",
    x = "Number of cells",
    y = "Spike fraction",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(out_dir, "spike_fraction_vs_ncells.png"),
  fraction_plot,
  width = 10,
  height = 5.8,
  dpi = 240
)

edge_plot <- ggplot(spike_df, aes(n_cells, mp_upper, color = dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  facet_wrap(~ fit_type, nrow = 1) +
  labs(
    title = "Fitted MP upper edge across cell counts",
    x = "Number of cells",
    y = "MP upper edge",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(out_dir, "mp_upper_edge_vs_ncells.png"),
  edge_plot,
  width = 10,
  height = 5.8,
  dpi = 240
)

compact_df <- spike_df[
  spike_df$fit_type == "Fixed gamma = p / n",
  c("dataset", "n_cells", "spike_count", "spike_count_for_plot", "spike_fraction_all", "mp_upper", "max_eigenvalue")
]
compact_df <- compact_df[order(compact_df$dataset, compact_df$n_cells), ]
print(compact_df, row.names = FALSE)
