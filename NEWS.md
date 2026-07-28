# scScale 0.0.3

- Added `scscale()` as the single user-facing entry point for the complete
  paired UMI-by-cell-number scaling workflow.
- Added `plot.scscale_model()` for UMI, cell-number, and joint-grid scaling
  plots.

# scScale 0.0.2

- `scscale_umi_mi()` now saves a long per-reference-spike `scaling` table.
- Added per-spike linear UMI parameters in `q_fit_rate`: intercept, slope,
  fitted R-squared, number of sampled rates, and number of detected rates.
- Added `q_X_rate_hat` to the long table and total-spike summaries in `q_total`
  and `q_total_rate_fit`.
