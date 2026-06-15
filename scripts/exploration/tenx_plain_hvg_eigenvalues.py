import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def read_tenx_matrix(data_dir):
    data_dir = Path(data_dir)
    matrix_file = data_dir / "matrix.mtx.gz"
    features_file = data_dir / "features.tsv.gz"
    barcodes_file = data_dir / "barcodes.tsv.gz"

    counts = scipy.io.mmread(matrix_file).tocsr()
    with gzip.open(features_file, "rt") as handle:
        features = pd.read_csv(handle, sep="\t", header=None)
    with gzip.open(barcodes_file, "rt") as handle:
        barcodes = [line.strip() for line in handle]
    return counts, features, barcodes


def normalize_log1p_gene_by_cell(counts, scale_factor=1e4):
    lib_size = np.asarray(counts.sum(axis=0)).ravel()
    lib_size[lib_size <= 0] = 1
    x = counts.tocsc(copy=True) @ sp.diags(scale_factor / lib_size)
    x = x.tocsr()
    x.data = np.log1p(x.data)
    return x


def sparse_row_mean_var(x):
    mean = np.asarray(x.mean(axis=1)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=1)).ravel()
    var = mean_sq - mean**2
    var[var < 0] = 0
    return mean, var


def covariance_eigenvalues(x_gene_by_cell):
    x_dense = x_gene_by_cell.toarray() if sp.issparse(x_gene_by_cell) else np.asarray(x_gene_by_cell)
    x_dense = x_dense - x_dense.mean(axis=1, keepdims=True)
    n_cells = x_dense.shape[1]
    singular_values = np.linalg.svd(x_dense, full_matrices=False, compute_uv=False)
    eigvals = (singular_values**2) / (n_cells - 1)
    if len(eigvals) < x_dense.shape[0]:
        eigvals = np.concatenate([eigvals, np.zeros(x_dense.shape[0] - len(eigvals))])
    eigvals = np.sort(eigvals)[::-1]
    eigvals[(eigvals < 0) & (eigvals > -1e-8)] = 0
    return eigvals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tenx_dir")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", required=True)
    parser.add_argument("--n-variable-genes", type=int, default=2000)
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--cell-counts", default=None)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts, features, barcodes = read_tenx_matrix(args.tenx_dir)
    rng = np.random.default_rng(args.seed)
    n_total_cells = counts.shape[1]
    cell_counts = None
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

    x = normalize_log1p_gene_by_cell(counts)
    _, gene_var = sparse_row_mean_var(x)
    n_genes = min(args.n_variable_genes, x.shape[0])
    gene_idx = np.sort(np.argsort(gene_var)[-n_genes:])
    x_hvg_full = x[gene_idx, :]

    eig_rows = []
    summary_rows = []
    for n_cells in cell_counts:
        if n_cells == n_total_cells:
            cell_idx = np.arange(n_total_cells)
        else:
            cell_idx = np.sort(rng.choice(n_total_cells, size=n_cells, replace=False))
        x_hvg = x_hvg_full[:, cell_idx]
        eigvals = covariance_eigenvalues(x_hvg)
        eig_rows.append(
            pd.DataFrame(
                {
                    "n_cells": x_hvg.shape[1],
                    "rank": np.arange(1, len(eigvals) + 1),
                    "eigenvalue": eigvals,
                    "log10_eigenvalue": np.log10(np.maximum(eigvals, np.finfo(float).tiny)),
                }
            )
        )
        summary_rows.append(
            {
                "dataset": args.dataset_name,
                "tenx_dir": args.tenx_dir,
                "n_total_cells": n_total_cells,
                "n_cells": x_hvg.shape[1],
                "n_genes_total": x.shape[0],
                "n_genes_used": x_hvg.shape[0],
                "hvg_method": "plain_log_normalized_variance",
                "max_eigenvalue": float(np.max(eigvals)),
                "median_eigenvalue": float(np.median(eigvals)),
                "mean_eigenvalue": float(np.mean(eigvals)),
            }
        )

    eig_df = pd.concat(eig_rows, ignore_index=True)
    eig_file = out_dir / f"{args.dataset_name}_plain_hvg_eigenvalues_by_cells.csv"
    eig_df.to_csv(eig_file, index=False)

    selected_features = features.iloc[gene_idx].copy()
    selected_features.to_csv(out_dir / f"{args.dataset_name}_plain_hvg_selected_features.csv", index=False)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / f"{args.dataset_name}_plain_hvg_eigenvalues_summary.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Wrote eigenvalues to {eig_file}")


if __name__ == "__main__":
    main()
