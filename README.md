# scScale

`scScale` is a small R package for single-cell scaling-law analysis. It fits
Gaussian spike models to count matrices and uses one low-rank subspace-alignment
mutual information formula to study matched RNA/ADT modalities.

The main workflow is intentionally small:

- `scscale_pair_fit()` fits the two modalities once and stores their aligned
  low-rank subspaces.
- `scscale_umi_mi()` refits the RNA spike model over UMI depth and computes the
  refitted theory curve plus the `I_infinity` bound.
- `scscale_cell_number_mi()` computes the cell-number scaling curve.
- `scscale_cell_number_by_umi_mi()` computes the joint cell-number by UMI
  scaling surface.
- `scscale_empirical_mi()` gives the empirical comparison curve.
- `scscale_low_rank_mi()` is the shared MI formula used by the theory helpers.

```r
library(scScale)

data(gse164378_3p_citeseq_hvg)
x <- gse164378_3p_citeseq_hvg$rna_counts
y <- gse164378_3p_citeseq_hvg$adt_counts

pair <- scscale_pair_fit(x, y)
umi <- scscale_umi_mi(pair, x, sampling_rates = c(0.25, 0.5, 1))
cells <- scscale_cell_number_mi(pair, n_grid = c(500, 1000, 2000, 4000))
joint <- scscale_cell_number_by_umi_mi(pair, umi, n_grid = c(500, 1000, 2000, 4000))
```

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

The rendered tutorial and reference files are available on GitHub Pages:

- [scScale tutorial index](https://chesleaz.github.io/scScale/)
- [Package manual PDF](https://chesleaz.github.io/scScale/scScale-manual.pdf)
- [GSE164378 RNA-ADT mutual information](https://chesleaz.github.io/scScale/tutorials/gse164378-rna-adt-mi.html)
- [GSE164378 RNA UMI scaling](https://chesleaz.github.io/scScale/tutorials/gse164378-umi-scaling.html)
- [GSE164378 Gaussian spike fit](https://chesleaz.github.io/scScale/tutorials/gse164378-scscale-fit.html)

The same examples are also installed with the package:

```r
library(scScale)
tutorial_dir <- system.file("tutorials", package = "scScale")
list.files(tutorial_dir, pattern = "\\.html$", full.names = TRUE)
```

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
