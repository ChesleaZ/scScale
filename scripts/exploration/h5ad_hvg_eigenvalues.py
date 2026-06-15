import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def sparse_mean_var(x):
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    var = mean_sq - mean**2
    var[var < 0] = 0
    return mean, var


def normalize_log1p(x, scale_factor=1e4):
    lib_size = np.asarray(x.sum(axis=1)).ravel()
    lib_size[lib_size <= 0] = 1
    if sp.issparse(x):
        x = x.tocsr(copy=True)
        scale = scale_factor / lib_size
        x = sp.diags(scale) @ x
        x.data = np.log1p(x.data)
        return x
    return np.log1p(x / lib_size[:, None] * scale_factor)


def dense_gene_cov_eigenvalues(x_cells_by_genes):
    if sp.issparse(x_cells_by_genes):
        x_genes_by_cells = x_cells_by_genes.T.tocsr()
        gene_means = np.asarray(x_genes_by_cells.mean(axis=1)).ravel()
        cov = (x_genes_by_cells @ x_genes_by_cells.T).toarray()
    else:
        x_genes_by_cells = np.asarray(x_cells_by_genes).T
        gene_means = x_genes_by_cells.mean(axis=1)
        cov = x_genes_by_cells @ x_genes_by_cells.T

    n_cells = x_cells_by_genes.shape[0]
    cov = (cov - n_cells * np.outer(gene_means, gene_means)) / (n_cells - 1)
    cov = (cov + cov.T) / 2
    eigvals = np.linalg.eigvalsh(cov)[::-1]
    eigvals[(eigvals < 0) & (eigvals > -1e-8)] = 0
    return eigvals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", default="h5ad_dataset")
    parser.add_argument("--n-variable-genes", type=int, default=2000)
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--cell-counts", default=None)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--assume-log-normalized", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    adata = ad.read_h5ad(args.h5ad)
    n_total_cells = adata.n_obs

    x = adata.X
    if not args.assume_log_normalized:
        x = normalize_log1p(x)

    _, gene_var = sparse_mean_var(x) if sp.issparse(x) else (x.mean(axis=0), x.var(axis=0))
    n_genes = min(args.n_variable_genes, x.shape[1])
    gene_idx = np.argsort(gene_var)[-n_genes:]
    gene_idx = np.sort(gene_idx)
    x_hvg_full = x[:, gene_idx]

    if args.cell_counts is not None:
        cell_counts = []
        for item in args.cell_counts.split(","):
            item = item.strip()
            if item.lower() in {"all", "inf"}:
                cell_counts.append(n_total_cells)
            elif item:
                cell_counts.append(min(int(item), n_total_cells))
        cell_counts = sorted(set(cell_counts))
    elif args.max_cells is not None and n_total_cells > args.max_cells:
        cell_counts = [args.max_cells]
    else:
        cell_counts = [n_total_cells]

    eig_rows = []
    summary_rows = []
    for n_cells in cell_counts:
        if n_cells == n_total_cells:
            cell_idx = np.arange(n_total_cells)
        else:
            cell_idx = np.sort(rng.choice(n_total_cells, size=n_cells, replace=False))
        x_hvg = x_hvg_full[cell_idx, :]
        eigvals = dense_gene_cov_eigenvalues(x_hvg)
        eig_rows.append(
            pd.DataFrame(
                {
                    "n_cells": x_hvg.shape[0],
                    "rank": np.arange(1, len(eigvals) + 1),
                    "eigenvalue": eigvals,
                    "log10_eigenvalue": np.log10(np.maximum(eigvals, np.finfo(float).tiny)),
                }
            )
        )
        summary_rows.append(
            {
                "dataset": args.dataset_name,
                "h5ad": args.h5ad,
                "n_total_cells": n_total_cells,
                "n_cells": x_hvg.shape[0],
                "n_genes_total": adata.n_vars,
                "n_genes_used": x_hvg.shape[1],
                "assume_log_normalized": args.assume_log_normalized,
                "max_eigenvalue": float(np.max(eigvals)),
                "median_eigenvalue": float(np.median(eigvals)),
                "mean_eigenvalue": float(np.mean(eigvals)),
            }
        )

    eig_df = pd.concat(eig_rows, ignore_index=True)
    eig_file = out_dir / f"{args.dataset_name}_hvg_eigenvalues_by_cells.csv"
    eig_df.to_csv(eig_file, index=False)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / f"{args.dataset_name}_hvg_eigenvalues_summary.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Wrote eigenvalues to {eig_file}")


if __name__ == "__main__":
    main()
