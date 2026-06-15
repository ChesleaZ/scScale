#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-rmt-sc")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/rmt-sc-cache")

import h5py
import matplotlib
import numpy as np
import pandas as pd
from scipy import optimize, sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = ROOT / (
    "reports/datasets_python_backup_20260609_124826/"
    "10x_PBMC_Multiome_3k/10x_PBMC_Multiome_3k_processed_pca.h5ad"
)
DEFAULT_OUT = ROOT / "outputs/exploration/pbmc_smallest_10x_scaling"


@dataclass(frozen=True)
class MPFit:
    gamma: float
    sigma2: float
    mp_lower: float
    mp_upper: float
    objective: float
    n_fit: int


def read_h5ad_csr_counts(path: Path, layer: str = "layers/counts") -> sparse.csr_matrix:
    with h5py.File(path, "r") as handle:
        group = handle[layer]
        shape = tuple(int(x) for x in group.attrs["shape"])
        counts = sparse.csr_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=shape,
        )
    return counts.astype(np.float64)


def log_normalize(counts: sparse.csr_matrix, scale_factor: float = 1e4) -> np.ndarray:
    lib_size = np.asarray(counts.sum(axis=1)).ravel()
    lib_size[lib_size <= 0] = 1.0
    normalized = counts.multiply(scale_factor / lib_size[:, None])
    return np.log1p(normalized.toarray())


def covariance_eigenvalues(cells_by_genes: np.ndarray) -> np.ndarray:
    x = cells_by_genes - cells_by_genes.mean(axis=0, keepdims=True)
    singular_values = np.linalg.svd(x, full_matrices=False, compute_uv=False)
    eigenvalues = singular_values**2 / max(cells_by_genes.shape[0] - 1, 1)
    return np.sort(eigenvalues[eigenvalues > 1e-12])[::-1]


def mp_positive_quantiles(gamma: float, sigma2: float, probs: np.ndarray, grid_n: int = 6000) -> np.ndarray:
    lower = sigma2 * (1.0 - math.sqrt(gamma)) ** 2
    upper = sigma2 * (1.0 + math.sqrt(gamma)) ** 2
    x = np.linspace(max(lower, 0.0) + 1e-12, upper - 1e-12, grid_n)
    density = np.zeros_like(x)
    inside = (x >= lower) & (x <= upper) & (x > 0)
    density[inside] = np.sqrt((upper - x[inside]) * (x[inside] - lower)) / (
        2.0 * math.pi * gamma * sigma2 * x[inside]
    )
    positive_mass = 1.0 / gamma if gamma > 1.0 else 1.0
    density = density / positive_mass
    dx = np.diff(x, append=x[-1] + (x[-1] - x[-2]))
    cdf = np.cumsum(density * dx)
    cdf = cdf / np.nanmax(cdf)
    return np.interp(probs, cdf, x)


def fit_free_gamma_mp(values: np.ndarray, gamma_start: float) -> MPFit:
    values = np.asarray(values)
    values = values[np.isfinite(values) & (values > 1e-10)]
    fit_probs = np.arange(0.05, 0.8001, 0.025)
    empirical_q = np.quantile(values, fit_probs)

    def objective(par: np.ndarray) -> float:
        gamma, sigma2 = np.exp(par)
        fitted_q = mp_positive_quantiles(gamma, sigma2, fit_probs)
        return float(np.mean((np.log(empirical_q) - np.log(fitted_q)) ** 2))

    starts = []
    for gamma in [gamma_start, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        for sigma2 in [np.median(values), np.mean(values), np.quantile(values, 0.75)]:
            starts.append((max(gamma, 1e-4), max(sigma2, np.finfo(float).eps)))

    bounds = [(math.log(1e-4), math.log(50.0)), (math.log(values.min() / 100.0), math.log(values.max() * 100.0))]
    best = None
    for start in starts:
        result = optimize.minimize(objective, np.log(start), method="L-BFGS-B", bounds=bounds, options={"maxiter": 800})
        if best is None or result.fun < best.fun:
            best = result

    gamma, sigma2 = np.exp(best.x)
    lower = sigma2 * (1.0 - math.sqrt(gamma)) ** 2
    upper = sigma2 * (1.0 + math.sqrt(gamma)) ** 2
    return MPFit(gamma=gamma, sigma2=sigma2, mp_lower=lower, mp_upper=upper, objective=float(best.fun), n_fit=len(values))


def fit_bulk_iterative(values: np.ndarray, gamma_theory: float, max_iter: int = 4) -> MPFit:
    positive = values[np.isfinite(values) & (values > 1e-10)]
    fit_values = positive
    fit = None
    for _ in range(max_iter):
        fit = fit_free_gamma_mp(fit_values, gamma_theory)
        next_values = positive[positive <= fit.mp_upper]
        if len(next_values) < 50 or len(next_values) == len(fit_values):
            break
        fit_values = next_values
    return MPFit(
        gamma=fit.gamma,
        sigma2=fit.sigma2,
        mp_lower=fit.mp_lower,
        mp_upper=fit.mp_upper,
        objective=fit.objective,
        n_fit=len(fit_values),
    )


def score_spectrum(values: np.ndarray, fit: MPFit, n_cells: int, n_genes: int) -> dict[str, float]:
    spikes = values[values > fit.mp_upper]
    excess = np.maximum(spikes / fit.mp_upper - 1.0, 0.0)
    top20 = values[: min(20, len(values))]
    total = values.sum()
    return {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "mp_gamma": fit.gamma,
        "mp_sigma2": fit.sigma2,
        "noise_edge": fit.mp_upper,
        "mp_objective": fit.objective,
        "n_fit": fit.n_fit,
        "n_spikes": len(spikes),
        "lambda1_over_edge": values[0] / fit.mp_upper,
        "top20_mean_over_edge": float(np.mean(top20 / fit.mp_upper)),
        "normalized_spike_excess": float(excess.sum()),
        "sqrt_normalized_spike_excess": float(math.sqrt(excess.sum())),
        "spike_excess_fraction": float(np.sum(spikes - fit.mp_upper) / total) if total > 0 else np.nan,
        "spike_variance_fraction": float(spikes.sum() / total) if total > 0 else np.nan,
        "total_variance": float(total),
    }


def analyze_matrix(counts: sparse.csr_matrix, rng: np.random.Generator, n_cells: int, depth_rate: float) -> tuple[dict, np.ndarray]:
    if n_cells < counts.shape[0]:
        selected = np.sort(rng.choice(counts.shape[0], size=n_cells, replace=False))
        sub_counts = counts[selected, :]
    else:
        sub_counts = counts

    if depth_rate < 1.0:
        coo = sub_counts.tocoo()
        thinned = rng.binomial(coo.data.astype(np.int64), depth_rate)
        keep = thinned > 0
        sub_counts = sparse.csr_matrix((thinned[keep], (coo.row[keep], coo.col[keep])), shape=coo.shape)

    log_norm = log_normalize(sub_counts)
    eig = covariance_eigenvalues(log_norm)
    fit = fit_bulk_iterative(eig, gamma_theory=sub_counts.shape[1] / sub_counts.shape[0])
    metrics = score_spectrum(eig, fit, sub_counts.shape[0], sub_counts.shape[1])
    umi = np.asarray(sub_counts.sum(axis=1)).ravel()
    metrics.update(
        {
            "depth_rate": depth_rate,
            "mean_umi_per_cell": float(np.mean(umi)),
            "median_umi_per_cell": float(np.median(umi)),
            "total_umi": float(np.sum(umi)),
        }
    )
    return metrics, eig


def fit_power_laws(metrics: pd.DataFrame, x_col: str, group_name: str) -> pd.DataFrame:
    metric_cols = [
        "noise_edge",
        "n_spikes",
        "lambda1_over_edge",
        "top20_mean_over_edge",
        "sqrt_normalized_spike_excess",
        "spike_excess_fraction",
        "spike_variance_fraction",
    ]
    rows = []
    for metric in metric_cols:
        d = metrics[[x_col, metric]].replace([np.inf, -np.inf], np.nan).dropna()
        d = d[(d[x_col] > 0) & (d[metric] > 0)]
        if len(d) < 3:
            continue
        x = np.log10(d[x_col].to_numpy())
        y = np.log10(d[metric].to_numpy())
        slope, intercept = np.polyfit(x, y, 1)
        pred = intercept + slope * x
        ss_res = np.sum((y - pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        rows.append(
            {
                "scaling_axis": group_name,
                "x_col": x_col,
                "metric": metric,
                "n_points": len(d),
                "x_min": float(d[x_col].min()),
                "x_max": float(d[x_col].max()),
                "log10_intercept": float(intercept),
                "exponent": float(slope),
                "r2_log10": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_scaling(metrics: pd.DataFrame, x_col: str, title: str, out_path: Path) -> None:
    plot_metrics = [
        "noise_edge",
        "n_spikes",
        "lambda1_over_edge",
        "top20_mean_over_edge",
        "sqrt_normalized_spike_excess",
        "spike_excess_fraction",
    ]
    fig, axes = plt.subplots(3, 2, figsize=(11, 10), constrained_layout=True)
    for ax, metric in zip(axes.ravel(), plot_metrics):
        d = metrics[[x_col, metric]].dropna()
        d = d[(d[x_col] > 0) & (d[metric] > 0)]
        ax.plot(d[x_col], d[metric], marker="o", linewidth=1.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(metric)
        ax.set_xlabel(x_col)
        ax.set_ylabel("value")
    fig.suptitle(title)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    out_dir = Path(os.environ.get("PBMC_SMALLEST_SCALING_OUT", DEFAULT_OUT))
    out_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = Path(os.environ.get("PBMC_SMALLEST_H5AD", DEFAULT_H5AD))
    seed = int(os.environ.get("PBMC_SMALLEST_SEED", "17"))
    rng = np.random.default_rng(seed)

    cell_counts = [int(x) for x in os.environ.get("PBMC_CELL_COUNTS", "300,500,800,1200,1800,2400,2699").split(",")]
    depth_rates = [float(x) for x in os.environ.get("PBMC_DEPTH_RATES", "0.1,0.2,0.35,0.5,0.75,1.0").split(",")]

    counts = read_h5ad_csr_counts(h5ad_path)
    cell_counts = sorted({min(c, counts.shape[0]) for c in cell_counts if c > 30})

    metric_rows = []
    eigen_rows = []

    for n_cells in cell_counts:
        print(f"cell scaling: n_cells={n_cells}")
        metrics, eig = analyze_matrix(counts, rng, n_cells=n_cells, depth_rate=1.0)
        metrics.update({"scaling_axis": "cell_number", "subsample_label": f"cells_{n_cells}"})
        metric_rows.append(metrics)
        eigen_rows.append(pd.DataFrame({"scaling_axis": "cell_number", "subsample_label": f"cells_{n_cells}", "rank": np.arange(1, len(eig) + 1), "eigenvalue": eig}))

    for rate in depth_rates:
        print(f"UMI scaling: depth_rate={rate}")
        metrics, eig = analyze_matrix(counts, rng, n_cells=counts.shape[0], depth_rate=rate)
        metrics.update({"scaling_axis": "umi_depth", "subsample_label": f"depth_{rate:g}"})
        metric_rows.append(metrics)
        eigen_rows.append(pd.DataFrame({"scaling_axis": "umi_depth", "subsample_label": f"depth_{rate:g}", "rank": np.arange(1, len(eig) + 1), "eigenvalue": eig}))

    metrics_df = pd.DataFrame(metric_rows)
    eig_df = pd.concat(eigen_rows, ignore_index=True)
    cell_metrics = metrics_df[metrics_df["scaling_axis"] == "cell_number"].copy()
    depth_metrics = metrics_df[metrics_df["scaling_axis"] == "umi_depth"].copy()

    fit_df = pd.concat(
        [
            fit_power_laws(cell_metrics, "n_cells", "cell_number"),
            fit_power_laws(depth_metrics, "mean_umi_per_cell", "umi_depth"),
        ],
        ignore_index=True,
    )

    metrics_df.to_csv(out_dir / "10x_pbmc_smallest_scaling_metrics.csv", index=False)
    eig_df.to_csv(out_dir / "10x_pbmc_smallest_scaling_eigenvalues.csv", index=False)
    fit_df.to_csv(out_dir / "10x_pbmc_smallest_scaling_fit_parameters.csv", index=False)
    plot_scaling(cell_metrics, "n_cells", "10x PBMC 3k: cell-number scaling", out_dir / "10x_pbmc_cell_number_scaling.png")
    plot_scaling(depth_metrics, "mean_umi_per_cell", "10x PBMC 3k: UMI-depth scaling", out_dir / "10x_pbmc_umi_depth_scaling.png")

    metadata = {
        "dataset_id": "10x_PBMC_Multiome_3k",
        "dataset_label": "10x PBMC Multiome 3k",
        "input_h5ad": str(h5ad_path.relative_to(ROOT) if h5ad_path.is_relative_to(ROOT) else h5ad_path),
        "counts_layer": "layers/counts",
        "n_cells_total": int(counts.shape[0]),
        "n_genes": int(counts.shape[1]),
        "seed": seed,
        "cell_counts": cell_counts,
        "depth_rates": depth_rates,
        "outputs": [
            "10x_pbmc_smallest_scaling_metrics.csv",
            "10x_pbmc_smallest_scaling_eigenvalues.csv",
            "10x_pbmc_smallest_scaling_fit_parameters.csv",
            "10x_pbmc_cell_number_scaling.png",
            "10x_pbmc_umi_depth_scaling.png",
        ],
    }
    (out_dir / "10x_pbmc_smallest_scaling_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
