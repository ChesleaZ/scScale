import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_zip_member_lines(zip_file, member):
    proc = subprocess.Popen(
        ["unzip", "-p", zip_file, member],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    if proc.stdout is None:
        raise RuntimeError("Could not open unzip stream")
    return proc


def normalize_log1p_csr(x, scale_factor=1e4):
    lib_size = np.asarray(x.sum(axis=1)).ravel()
    lib_size[lib_size <= 0] = 1
    x = x.tocsr(copy=True)
    x = sp.diags(scale_factor / lib_size) @ x
    x.data = np.log1p(x.data)
    return x


def sparse_col_mean_var(x):
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    var = mean_sq - mean**2
    var[var < 0] = 0
    return mean, var


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_file")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", default="parse_1m_pbmc_t1d")
    parser.add_argument("--sample-cells", type=int, default=100000)
    parser.add_argument("--n-variable-genes", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    proc = read_zip_member_lines(args.zip_file, "DGE_1M_PBMC.mtx")
    header = proc.stdout.readline().strip()
    if not header.startswith("%%MatrixMarket"):
        raise ValueError(f"Unexpected MatrixMarket header: {header}")

    line = proc.stdout.readline().strip()
    while line.startswith("%"):
        line = proc.stdout.readline().strip()
    n_cells, n_genes, n_nz = map(int, line.split())

    sample_n = min(args.sample_cells, n_cells)
    sampled_cells = np.sort(rng.choice(n_cells, size=sample_n, replace=False) + 1)
    sampled_set = set(sampled_cells.tolist())
    sampled_map = {cell: i for i, cell in enumerate(sampled_cells.tolist())}

    rows = []
    cols = []
    data = []
    kept_nz = 0
    for line in proc.stdout:
      if not line:
          continue
      i, j, v = line.split()
      i = int(i)
      if i not in sampled_set:
          continue
      rows.append(sampled_map[i])
      cols.append(int(j) - 1)
      data.append(float(v))
      kept_nz += 1

    status = proc.wait()
    if status != 0:
        raise RuntimeError(f"unzip exited with status {status}")

    counts = sp.coo_matrix((data, (rows, cols)), shape=(sample_n, n_genes)).tocsr()
    counts.eliminate_zeros()

    x = normalize_log1p_csr(counts)
    _, gene_var = sparse_col_mean_var(x)
    gene_idx = np.sort(np.argsort(gene_var)[-args.n_variable_genes:])
    x_hvg = x[:, gene_idx]

    x_genes_by_cells = x_hvg.T.tocsr()
    gene_means = np.asarray(x_genes_by_cells.mean(axis=1)).ravel()
    cov = (x_genes_by_cells @ x_genes_by_cells.T).toarray()
    cov = (cov - sample_n * np.outer(gene_means, gene_means)) / (sample_n - 1)
    cov = (cov + cov.T) / 2
    eigvals = np.linalg.eigvalsh(cov)[::-1]
    eigvals[(eigvals < 0) & (eigvals > -1e-8)] = 0

    eig_df = pd.DataFrame(
        {
            "n_cells": sample_n,
            "rank": np.arange(1, len(eigvals) + 1),
            "eigenvalue": eigvals,
            "log10_eigenvalue": np.log10(np.maximum(eigvals, np.finfo(float).tiny)),
        }
    )
    eig_file = out_dir / f"{args.dataset_name}_hvg_eigenvalues_by_cells.csv"
    eig_df.to_csv(eig_file, index=False)

    meta_proc = read_zip_member_lines(args.zip_file, "cell_metadata_1M_PBMC.csv")
    metadata = pd.read_csv(meta_proc.stdout, usecols=["gene_count", "tscp_count", "read_count"])
    meta_status = meta_proc.wait()
    if meta_status != 0:
        raise RuntimeError(f"unzip metadata stream exited with status {meta_status}")
    sampled_meta = metadata.iloc[sampled_cells - 1]
    summary = pd.DataFrame(
        {
            "dataset": [args.dataset_name],
            "zip_file": [args.zip_file],
            "n_total_cells": [n_cells],
            "n_total_genes": [n_genes],
            "n_total_nnz": [n_nz],
            "sampled_cells": [sample_n],
            "sampled_nnz": [kept_nz],
            "n_genes_used": [len(gene_idx)],
            "median_gene_count": [float(sampled_meta["gene_count"].median())],
            "median_tscp_count": [float(sampled_meta["tscp_count"].median())],
            "mean_tscp_count": [float(sampled_meta["tscp_count"].mean())],
            "max_eigenvalue": [float(np.max(eigvals))],
            "median_eigenvalue": [float(np.median(eigvals))],
            "mean_eigenvalue": [float(np.mean(eigvals))],
        }
    )
    summary_file = out_dir / f"{args.dataset_name}_hvg_eigenvalues_summary.csv"
    summary.to_csv(summary_file, index=False)

    print(summary.to_string(index=False))
    print(f"Wrote eigenvalues to {eig_file}")


if __name__ == "__main__":
    main()
