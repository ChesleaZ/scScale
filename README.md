# scScale

Small R package for object-style single-cell scaling-law analysis.

The public workflow is intentionally narrow:

1. `scscale_fit()` fits the Gaussian spike model for one feature-by-cell matrix.
2. `scscale_mi()` uses a second fitted object as the target and computes
   spectral mutual information, optionally across cell-number and UMI grids.

```r
library(scScale)
data(gse164378_3p_citeseq_hvg)

rna <- gse164378_3p_citeseq_hvg$rna_counts
adt <- gse164378_3p_citeseq_hvg$adt_counts

rna_fit <- scscale_fit(
  rna,
  r = 12,
  fit_umi = TRUE,
  sampling_rates = c(0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00)
)

adt_fit <- scscale_fit(adt, r = 12, fit_umi = FALSE)

mi <- scscale_mi(
  rna_fit,
  adt_fit,
  n_grid = c(500, 1000, 2000, 3000, 4000),
  sampling_rates = c(0.35, 0.50, 0.70, 0.85, 1.00)
)

mi$I_theory
mi$I_infinity
mi$grid
```

`scscale_fit()` returns the fitted spectrum, MP bulk fit, spike table, finite
recoverability `theta_X`, infinite-cell recoverability `theta_infinity`, and
optionally UMI scaling. The spike table contains the fitted `d2_X`, `q_X`, MP
edge, and related parameters.

The package examples are:

- `gse164378_3p_citeseq_hvg`: 4,000 matched 3' CITE-seq PBMC cells with 2,000
  RNA HVGs and 228 ADT features.
- `gse123025_myeloid_hvg`: 2,000 HVGs by 1,922 myeloid cells for a compact RNA
  spike-model example.

Tutorials are available under:

```r
system.file("tutorials", package = "scScale")
```

The older eigenvalue-first API is kept locally under `local/legacy_scScale/`
for provenance, but is no longer part of the GitHub package surface.
