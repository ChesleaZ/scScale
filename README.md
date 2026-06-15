# scScale

`scScale` is an R package for single-cell scaling-law analyses. It includes
median Marchenko-Pastur noise calibration, BBP spike inversion, spectral mutual
information utilities, and batch-effect scaling-law fits.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("ChesleaZ/Scaling-law-summer-version-")
```

You can also install with `pak`:

```r
install.packages("pak")
pak::pak("ChesleaZ/Scaling-law-summer-version-")
```

## Quick Start

```r
library(scScale)

ref_ev <- find_eigenvalues(ref_counts, n_features = 2000)
cur_ev <- find_eigenvalues(cur_counts, n_features = 2000)

ref_noise <- fit_noise(ref_ev)
ref <- fit_spikes(ref_ev, noise = ref_noise, r = 5)
cur <- fit_spikes(cur_ev, noise = ref_noise, r = 5)

mi_theory(cur, ref, side = "cells", r = 5)$mi
```

If eigenvalues are already precomputed:

```r
fit <- mi_theory_from_eigenvalues(
  cur_eigenvalues = cur_mu,
  ref_eigenvalues = ref_mu,
  n_cur = 600,
  p_cur = 300,
  side = "cells",
  r = 5
)
```

Batch-effect scaling fits can be run from replicate summaries:

```r
summary_df <- data.frame(
  m_batch = c(2, 3, 5, 8, 12),
  mean_I_bio_norm = c(0.12, 0.15, 0.18, 0.20, 0.215)
)

fit_batch_effect_scaling(summary_df, law = "batch_number", min_points = 4)
```

## Command Line

From precomputed eigenvalues:

```bash
Rscript inst/scripts/run_from_eigenvalues.R \
  --cur-eigenvalues=cur.csv \
  --ref-eigenvalues=ref.csv \
  --n-cur=600 \
  --p-cur=300 \
  --side=cells \
  --out-dir=outputs/scscale_spectral_mi
```

From count CSVs:

```bash
Rscript inst/scripts/run_from_eigenvalues.R \
  --cur-counts=cur_counts.csv \
  --ref-counts=ref_counts.csv \
  --n-features=2000 \
  --transform=pearson \
  --side=cells \
  --out-dir=outputs/scscale_spectral_mi
```

## Tutorials

Rendered tutorials are in [`docs/`](docs/):

- [Tutorial index](docs/scScale-tutorial.html)
- [Cell-number scaling law](docs/cell-number-scaling-law.html)
- [UMI scaling law](docs/umi-scaling-law.html)
- [Batch-number scaling law](docs/batch-number-scaling-law.html)

The GitHub Pages website entry point is [`docs/index.html`](docs/index.html).
When Pages is enabled for this repository using the `main` branch and `/docs`
folder, the tutorials will be available at:

<https://chesleaz.github.io/Scaling-law-summer-version-/>
