\name{scscale_representation}
\alias{scscale_representation}
\alias{scscale_pca_representation}
\alias{scscale_representation_umi_fit}
\alias{scscale_representation_cell_fit}
\title{Scaling Laws for PCA and External Representations}
\description{
Constructs a common representation object for PCA or externally generated
cell embeddings, fits the UMI-anchored effective model, and optionally
calibrates its held-out cell-number law.
}
\usage{
scscale_representation(
  reference,
  umi = NULL,
  umi_values = NULL,
  name = "external",
  rank = ncol(reference),
  rmt_dimension = ncol(reference)
)

scscale_pca_representation(
  reference_counts,
  umi_counts = NULL,
  umi_values = NULL,
  n_features = 2000,
  rank = 10,
  min_cells = 10,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE
)

scscale_representation_umi_fit(
  representation,
  target,
  target_rank = 12,
  target_input = c("counts", "normalized", "representation"),
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  weighting = c("inverse_variance", "equal"),
  n_starts = 1L,
  seed = 1,
  I_infinity_max_factor = 2,
  control = list(maxit = 5000)
)

scscale_representation_cell_fit(
  object,
  cell_numbers,
  replicates = 50,
  evaluation_cells = 5000,
  seed = 1,
  weighting = c("inverse_variance", "equal"),
  relative_floor = 1e-8
)
}
\arguments{
  \item{reference}{Full-depth cell-by-dimension external representation. Cell
  IDs must be unique row names.}
  \item{umi}{List of UMI-level external representations in the same format.}
  \item{umi_values}{Absolute median UMI per cell for each list element. Numeric
  list names are used when this is omitted. Repeated values represent
  replicates.}
  \item{name}{Display name for an external representation.}
  \item{rank}{Number of leading representation dimensions retained.}
  \item{rmt_dimension}{Dimension \eqn{d} used by the effective scaling model.
  This is separate from the retained rank.}
  \item{reference_counts}{Full-depth feature-by-cell count matrix.}
  \item{umi_counts}{List of UMI-thinned feature-by-cell count matrices.}
  \item{n_features}{Number of HVGs selected once from the reference counts and
  reused at every UMI level.}
  \item{min_cells}{Minimum detected cells used during HVG selection.}
  \item{target_depth,count_transform,center,scale}{Count normalization options
  shared with the standard scScale preprocessing path.}
  \item{representation}{A \code{scscale_representation} object.}
  \item{target}{A feature-by-cell or cell-by-feature target matrix, or a named
  numeric or categorical cell-level target.}
  \item{target_rank}{Number of target subspace dimensions retained.}
  \item{target_input}{Whether a matrix target contains counts, normalized
  values, or an existing representation. Ignored for vector targets.}
  \item{weighting}{Weighting passed to the effective UMI or cell calibration.}
  \item{n_starts,seed,I_infinity_max_factor,control}{Options passed to
  \code{scscale_effective_fit()}. The seed also defines the held-out split.}
  \item{object}{A fitted \code{scscale_representation_umi_fit} object.}
  \item{cell_numbers}{Training-set sizes used to learn paired CCA directions.}
  \item{replicates}{Number of nested random training permutations.}
  \item{evaluation_cells}{Number of cells reserved once for evaluating every
  cell-number level.}
  \item{relative_floor}{Relative covariance eigenvalue floor used only while
  whitening the training data for CCA.}
}
\details{
Both constructors return the same object. The PCA constructor is only an
adapter: it selects HVGs from the reference matrix, reuses them at all UMI
levels, applies standard scScale normalization, and supplies the resulting
cell subspaces to the common downstream code.

For UMI scaling, the target subspace is learned once from the full target and
restricted to the cells present at each UMI level. Each representation is
truncated independently to \code{rank}; empirical subspace MI is then passed
unchanged to \code{scscale_effective_fit()}.

For cell-number scaling, one fixed evaluation set is held out. At every cell
number, paired CCA directions are learned from a nested training subsample and
MI is evaluated on the same held-out cells. The UMI coefficients remain fixed;
\code{scscale_effective_cell_calibrate_B()} estimates only the cell-number
scale.
}
\value{
The constructors return a \code{scscale_representation}. The UMI fitting
function returns a \code{scscale_representation_umi_fit} that also inherits
from \code{scscale_effective_fit}. The cell function adds empirical replicate
points, component correlations, split metadata, and cell calibration.
}
\examples{
set.seed(1)
cells <- paste0("c", 1:60)
embedding <- matrix(rnorm(60 * 8), 60, 8, dimnames = list(cells, NULL))
umi <- lapply(1:5, function(i) embedding + matrix(rnorm(60 * 8, sd = 1 / i), 60, 8))
umi <- lapply(umi, function(x) { rownames(x) <- cells; x })
rep <- scscale_representation(
  embedding, umi, umi_values = c(100, 200, 400, 800, 1600), rank = 5
)
labels <- setNames(rep(c("A", "B"), 30), cells)
fit <- scscale_representation_umi_fit(rep, labels, target_rank = 1)
}
