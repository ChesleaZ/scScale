# scScale

`scScale` is a small R package for single-cell scaling-law analysis. It
supports a calibrated Gaussian spike model and a separate direct empirical
curve fit for matched RNA/ADT modalities.

The canonical spike workflow separates calibration from evaluation:

- `scscale_spike_fit()` performs all data-dependent work once: X and Y
  spike fitting, rank selection, observed subspace alignment, reference
  deattenuation, fixed-reference projections at each UMI depth, and per-spike
  zero-anchored ordinary least-squares \(q(U)=bU\) fitting.
- `scscale_cell_scaling()` evaluates \(I(n,1)\) from the saved fit.
- `scscale_umi_scaling()` evaluates \(I(n_{\mathrm{ref}},\rho)\).
- `scscale_joint_scaling()` evaluates the joint \(I(n,\rho)\) grid.

All three scaling functions use the same stored \(q(U)\) depth trajectories
and alignment parameters. Cell number enters only through
\(c_X(n)=p_X/n\); the functions then compute \(\theta_X(n,U)\) and MI.

```r
library(scScale)

data(gse164378_3p_citeseq_hvg)
x <- gse164378_3p_citeseq_hvg$rna_counts
y <- gse164378_3p_citeseq_hvg$adt_counts

fit <- scscale_spike_fit(x, y)

fit$x$parameters
fit$y$parameters
fit$alignment$P
fit$alignment$K

cell <- scscale_cell_scaling(fit, n = c(500, 1000, 2000, 4000))
umi <- scscale_umi_scaling(fit, sampling_rate = seq(0.1, 1, length.out = 50))
joint <- scscale_joint_scaling(
  fit,
  n = c(500, 1000, 2000, 4000),
  sampling_rate = c(0.2, 0.5, 1)
)
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
