# RMT_sc

Random matrix theory experiments for single-cell covariance spectra.

## Current Layout

- `data/Jurkat/`: Jurkat 10x matrix input.
- `data/K562/`: K562 10x matrix input.
- `scripts/exploration/`: one simple bucket for exploratory scripts, including the MP analysis pipeline work.
- `outputs/exploration/`: exploratory CSVs and plots, kept for reference.

## Current Pipeline Direction

The active analysis is exploratory MP work, moving toward mixture models as the core method:

1. detect spike eigenvalues,
2. estimate residual noise variance from non-spike eigenvalues,
3. compute the MP law using known aspect ratio and residual variance,
4. test whether one exact MP or a mixture of exact MPs explains the residual bulk.

Start with:

```sh
Rscript scripts/exploration/simulated_spiked_mp_pipeline.R
```
