# Scripts

This folder is intentionally simple during exploration.

## exploration

- All exploratory scripts live in `scripts/exploration/`.
- The MP analysis pipeline from the current exploration lives here too.
- These scripts are kept for reference, but most are not expected to become reusable production code.
- `fit_two_mp_mixture_em.R` is the first-pass EM prototype for two MP components plus an upper-tail spike/outlier component.
- `fit_two_mp_free_lambda_em.R` is the less constrained EM prototype where the two MP lambda parameters are fitted directly.
- `fit_gaussian_init_mp_em.R` restarts the mixture fit by first fitting a two-Gaussian-plus-spike model, then using that split to initialize MP refinement.
- `fit_gmm_spectrum.R` fits the current two-mode GMM plus fixed upper-tail outlier baseline.
- `fit_mp_from_gmm_init_em.R` removes the GMM outlier component, initializes two MP components from the two GMM modes, then refines with EM.
- `fit_shifted_mp_from_gmm_init_em.R` adds a fitted location shift (`delta`) to each MP component after GMM initialization and outlier removal.
- `fit_shifted_mp_mixture_from_gmm_em.R` is the current shifted-MP mixture EM attempt using the validated single-component M-step with `lambda <= 1`.
- `plot_shifted_mp_gmm_initialization.R` shows the shifted-MP distribution implied by fixed GMM responsibilities before any EM updates.
