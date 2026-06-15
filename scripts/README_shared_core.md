# Shared Scaling-Law Scripts

This folder separates reusable code from one-off exploration scripts.

## Core Modules

- `scripts/core/empirical_spectrum.R`
  - covariance eigenvalue calculation from cell-by-feature matrices
  - empirical CDF helpers
  - eigenvalue summary helpers

- `scripts/core/mp_bbp.R`
  - Marchenko-Pastur density, edges, quantiles
  - median fixed-gamma MP bulk fitting
  - fixed/free-gamma quantile MP bulk fitting for comparisons
  - iterative spike removal
  - canonical median MP fit: `fit_mp_median()`
  - BBP forward and inverse spike maps
  - combined `fit_bbp_spikes()` with `spike_strength_d`

- `scripts/core/spectral_overlap_mi.R`
  - empirical cell/gene subspace-overlap MI
  - theoretical overlap MI from BBP spike strengths

Source these files from other scripts:

```r
source("scripts/core/empirical_spectrum.R")
source("scripts/core/mp_bbp.R")
source("scripts/core/spectral_overlap_mi.R")
```

## Utility CLIs

Compute empirical covariance eigenvalues from a CSV/TSV matrix:

```sh
Rscript scripts/utils/compute_empirical_spectrum.R \
  --input=data/example_matrix.csv \
  --input-type=matrix \
  --out=outputs/example/eigenvalues.csv \
  --orientation=cells_by_features
```

Use precomputed PCA scores:

```sh
Rscript scripts/utils/compute_empirical_spectrum.R \
  --input=outputs/example/pca_scores.csv \
  --input-type=pca_scores \
  --out=outputs/example/eigenvalues.csv
```

Use PCA already stored in a Seurat object:

```sh
Rscript scripts/utils/compute_empirical_spectrum.R \
  --input=outputs/example/seurat_processed.rds \
  --input-type=seurat \
  --reduction=pca \
  --out=outputs/example/eigenvalues.csv
```

Pass through already-computed eigenvalues:

```sh
Rscript scripts/utils/compute_empirical_spectrum.R \
  --input=outputs/example/existing_eigenvalues.csv \
  --input-type=eigenvalues \
  --out=outputs/example/eigenvalues.csv
```

Fit MP bulk and BBP spikes from an eigenvalue CSV:

```sh
Rscript scripts/utils/run_bbp_fit_from_eigenvalues.R \
  --eigenvalues=outputs/example/eigenvalues.csv \
  --out-dir=outputs/example/bbp_fit \
  --dataset=example
```

The utility writes:

- `bbp_bulk_fit.csv`
- `bbp_spikes.csv`
- `mp_bulk_curve.csv`

For direct use inside analysis scripts, prefer:

```r
mp_fit <- fit_mp_median(
  eigenvalues,
  n_cells = n_cells,
  n_features = n_features,
  edge_buffer = 1.01
)
```

For theoretical spectral-overlap MI:

```r
cur <- fit_bbp_spikes(cur_eigenvalues, n_cells = n, n_features = p)
ref <- fit_bbp_spikes(ref_eigenvalues, n_cells = n, n_features = p)
mi <- theoretical_overlap_mi_from_spike_tables(cur$spikes, ref$spikes, side = "cells", r = 5)
```

## Notes

The GitHub repo linked in the notes was not accessible from this session without credentials, so these modules are factored from the local exploration scripts and rewritten into project-specific reusable functions.

For MP/BBP fitting, a full or nearly full PCA spectrum is much better than a
truncated 30-PC/50-PC embedding. The eigenvalue output includes
`spectrum_is_truncated` when the script can infer it. A truncated PCA is still
usable for quick diagnostics, but the fitted MP bulk can be biased because most
noise eigenvalues are missing.
