# scScale

Small R package-style module for single-cell scaling-law analyses.

The implementation includes median-MP/BBP spectral mutual information
utilities, plus batch-effect scaling-law fitting adapted from the BMMC CITE-seq
validation workflow in `ChesleaZ/Scaling-law-summer-version-`.

```r
library(scScale, lib.loc = ".rlib")

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

From the command line:

```bash
Rscript packages/scScale/inst/scripts/run_from_eigenvalues.R \
  --cur-eigenvalues=cur.csv \
  --ref-eigenvalues=ref.csv \
  --n-cur=600 \
  --p-cur=300 \
  --side=cells \
  --out-dir=outputs/scscale_spectral_mi
```

Or from count CSVs directly:

```bash
Rscript packages/scScale/inst/scripts/run_from_eigenvalues.R \
  --cur-counts=cur_counts.csv \
  --ref-counts=ref_counts.csv \
  --n-features=2000 \
  --transform=pearson \
  --side=cells \
  --out-dir=outputs/scscale_spectral_mi
```

Batch-effect scaling fits can be run from replicate summaries:

```r
summary_df <- data.frame(
  m_batch = c(2, 3, 5, 8, 12),
  mean_I_bio_norm = c(0.12, 0.15, 0.18, 0.20, 0.215)
)

fit_batch_effect_scaling(summary_df, law = "batch_number", min_points = 4)
```

Rendered tutorials are available in `packages/scScale/docs/`:

- `scScale-tutorial.html`: tutorial index
- `cell-number-scaling-law.html`
- `umi-scaling-law.html`
- `batch-number-scaling-law.html`

The tutorials start from a compact real Jurkat 10x count block derived from
`data/Jurkat/sample_filtered_feature_bc_matrix`, then vary cells, UMI depth, or
batch count.

The high-level API returns ordinary R-style fit objects:

```r
fit <- fit_cell_scaling(counts, cell_grid = c(500, 1000, 2000, 5000))
coef(fit)
predict(fit, data.frame(n_cells = c(10000, 50000)))
plot(fit)
```
