# scScale

`scScale` is a small R package for object-style single-cell scaling-law
analysis. It fits Gaussian spike models to count matrices and uses the fitted
spike parameters to study sequencing-depth scaling, cell-number scaling, and
spectral mutual information between modalities.

The package surface is intentionally narrow:

- `scscale_fit()` fits the Gaussian spike model for one feature-by-cell matrix.
- `scscale_mi()` compares two fitted objects and computes spectral mutual
  information across optional cell-number and UMI grids.

## Installation

Install the GitHub version with `remotes`:

```r
install.packages("remotes")
remotes::install_github("ChesleaZ/scScale")
```

or with `pak`:

```r
install.packages("pak")
pak::pak("ChesleaZ/scScale")
```

## Tutorials

The worked examples live in the installed package tutorials:

```r
library(scScale)
tutorial_dir <- system.file("tutorials", package = "scScale")
list.files(tutorial_dir, pattern = "\\.html$", full.names = TRUE)
```

Open the tutorial HTML files for the full analysis workflow, including Gaussian
spike fitting, UMI scaling, and CITE-seq mutual information.

## Example Data

- `gse164378_3p_citeseq_hvg`: 4,000 matched 3' CITE-seq PBMC cells with 2,000
  RNA HVGs and 228 ADT features.
- `gse123025_myeloid_hvg`: 2,000 HVGs by 1,922 myeloid cells for a compact RNA
  spike-model example.

## Documentation

```r
system.file("doc/scScale-manual.pdf", package = "scScale")
```

The older eigenvalue-first API is kept locally under `local/legacy_scScale/`
for provenance, but is no longer part of the GitHub package surface.
