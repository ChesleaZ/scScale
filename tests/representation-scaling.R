library(scScale)

set.seed(41)
n <- 90L
p <- 30L
cells <- paste0("cell_", seq_len(n))
features <- paste0("gene_", seq_len(p))
latent <- matrix(rnorm(3L * n), 3L, n)
loadings <- matrix(rnorm(p * 3L, sd = 0.35), p, 3L)
rate <- exp(1.2 + loadings %*% latent)
counts <- matrix(rpois(length(rate), rate), p, n,
  dimnames = list(features, cells))
target_rate <- exp(1 + matrix(rnorm(6L * 3L, sd = 0.4), 6L, 3L) %*% latent)
target <- matrix(rpois(length(target_rate), target_rate), 6L, n,
  dimnames = list(paste0("adt_", 1:6), cells))
fractions <- c(0.08, 0.15, 0.3, 0.5, 0.75, 1)
umi_counts <- lapply(seq_along(fractions), function(i) {
  scscale_downsample_counts_fraction(counts, fractions[[i]], seed = 100 + i)
})
umi_values <- fractions * 2000

pca <- scscale_pca_representation(
  counts,
  umi_counts = umi_counts,
  umi_values = umi_values,
  n_features = 20,
  rank = 5,
  min_cells = 1
)
external <- scscale_representation(
  reference = pca$reference,
  umi = pca$umi,
  umi_values = pca$umi_values,
  name = "PCA supplied externally",
  rank = pca$rank,
  rmt_dimension = pca$rmt_dimension
)
stopifnot(
  inherits(pca, "scscale_representation"),
  identical(pca$source, "PCA"),
  length(pca$features) == 20L,
  ncol(pca$reference) == 5L,
  length(pca$umi) == length(fractions)
)

pca_fit <- suppressWarnings(scscale_representation_umi_fit(
  pca, target, target_rank = 3, n_starts = 1, seed = 7
))
external_fit <- suppressWarnings(scscale_representation_umi_fit(
  external, target, target_rank = 3, n_starts = 1, seed = 7
))
stopifnot(
  inherits(pca_fit, "scscale_representation_umi_fit"),
  isTRUE(all.equal(
    pca_fit$umi_points$I_empirical,
    external_fit$umi_points$I_empirical,
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    unname(pca_fit$coefficients),
    unname(external_fit$coefficients),
    tolerance = 1e-8
  ))
)

cell_fit <- scscale_representation_cell_fit(
  pca_fit,
  cell_numbers = c(10, 25, 50),
  replicates = 2,
  evaluation_cells = 25,
  seed = 19,
  weighting = "equal"
)
stopifnot(
  inherits(cell_fit, "scscale_representation_cell_fit"),
  "B_reference" %in% names(cell_fit$coefficients),
  nrow(cell_fit$cell_points) == 2L * 4L,
  length(cell_fit$cell_design$evaluation_cells) == 25L,
  !length(intersect(
    cell_fit$cell_design$evaluation_cells,
    cell_fit$cell_design$training_pool
  )),
  all(is.finite(cell_fit$cell_points$I_empirical))
)
