library(scScale)

set.seed(42)
cells <- paste0("cell_", seq_len(60))
X <- matrix(rnorm(8 * 60), nrow = 8, dimnames = list(NULL, cells))
Y <- matrix(rnorm(7 * 60), nrow = 7, dimnames = list(NULL, cells))

mi <- scscale_empirical_mi(
  X,
  Y,
  input = "normalized",
  target_input = "normalized",
  x_rank = 3,
  y_rank = 5,
  use_irlba = FALSE,
  store_subspaces = TRUE
)
stopifnot(
  mi$x_rank == 3L,
  mi$y_rank == 5L,
  ncol(mi$z_X) == 3L,
  ncol(mi$z_Y) == 5L
)

grid <- scscale_empirical_inu_grid(
  X,
  Y,
  n_grid = c(30, 60),
  sampling_rates = 1,
  n_replicates = 1,
  input = "normalized",
  target_input = "normalized",
  x_rank = 3,
  y_rank = 5,
  use_irlba = FALSE
)
stopifnot(
  all(grid$x_rank == 3L),
  all(grid$y_rank == 5L)
)

# Existing shared-rank calls remain unchanged.
shared <- scscale_empirical_mi(
  X,
  Y,
  r = 4,
  input = "normalized",
  target_input = "normalized",
  use_irlba = FALSE
)
stopifnot(shared$x_rank == 4L, shared$y_rank == 4L)
