import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def open_zip_member(zip_file, member):
    proc = subprocess.Popen(
        ["unzip", "-p", zip_file, member],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    if proc.stdout is None:
        raise RuntimeError(f"Could not stream {member} from {zip_file}")
    return proc


def normalize_log1p_cells_by_genes(counts, scale_factor=1e4):
    lib_size = np.asarray(counts.sum(axis=1)).ravel()
    lib_size[lib_size <= 0] = 1
    x = counts.tocsr(copy=True)
    x = sp.diags(scale_factor / lib_size) @ x
    x.data = np.log1p(x.data)
    return x


def sparse_col_var(x):
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    var = mean_sq - mean**2
    var[var < 0] = 0
    return var


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_file")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", default="parse_pbmc")
    parser.add_argument("--sample-cells", type=int, default=5000)
    parser.add_argument("--n-variable-genes", type=int, default=400)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    proc = open_zip_member(args.zip_file, "DGE_1M_PBMC.mtx")
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
        cell, gene, value = line.split()
        cell = int(cell)
        if cell not in sampled_set:
            continue
        rows.append(sampled_map[cell])
        cols.append(int(gene) - 1)
        data.append(int(float(value)))
        kept_nz += 1

    status = proc.wait()
    if status != 0:
        raise RuntimeError(f"unzip DGE stream exited with status {status}")

    counts = sp.coo_matrix(
        (data, (rows, cols)),
        shape=(sample_n, n_genes),
        dtype=np.int32,
    ).tocsr()
    counts.eliminate_zeros()

    x = normalize_log1p_cells_by_genes(counts)
    gene_var = sparse_col_var(x)
    n_hvg = min(args.n_variable_genes, n_genes)
    gene_idx = np.sort(np.argsort(gene_var)[-n_hvg:])
    counts_hvg = counts[:, gene_idx].T.toarray()

    gene_proc = open_zip_member(args.zip_file, "all_genes_1M_PBMC.csv")
    genes = pd.read_csv(gene_proc.stdout)
    gene_status = gene_proc.wait()
    if gene_status != 0:
        raise RuntimeError(f"unzip gene stream exited with status {gene_status}")

    meta_proc = open_zip_member(args.zip_file, "cell_metadata_1M_PBMC.csv")
    metadata = pd.read_csv(meta_proc.stdout)
    meta_status = meta_proc.wait()
    if meta_status != 0:
        raise RuntimeError(f"unzip metadata stream exited with status {meta_status}")

    selected_features = genes.iloc[gene_idx].copy()
    selected_features.insert(0, "matrix_col", gene_idx + 1)
    selected_features.to_csv(out_dir / f"{args.dataset_name}_glmpca_selected_features.csv", index=False)

    selected_cells = metadata.iloc[sampled_cells - 1].copy()
    selected_cells.insert(0, "matrix_row", sampled_cells)
    selected_cells.to_csv(out_dir / f"{args.dataset_name}_glmpca_selected_cells.csv", index=False)

    gene_labels = selected_features["gene_name"].fillna(selected_features["gene_id"]).astype(str)
    gene_labels = [f"{label}_{idx}" for label, idx in zip(gene_labels, gene_idx + 1)]
    cell_labels = [f"cell_{i}" for i in sampled_cells]
    count_df = pd.DataFrame(counts_hvg, index=gene_labels, columns=cell_labels)
    count_file = out_dir / f"{args.dataset_name}_glmpca_counts_hvg.csv"
    count_df.to_csv(count_file)

    summary = pd.DataFrame(
        {
            "dataset": [args.dataset_name],
            "zip_file": [args.zip_file],
            "n_total_cells": [n_cells],
            "n_total_genes": [n_genes],
            "n_total_nnz": [n_nz],
            "sampled_cells": [sample_n],
            "sampled_nnz": [kept_nz],
            "n_genes_used": [n_hvg],
            "mean_umi_per_cell": [float(np.asarray(counts.sum(axis=1)).mean())],
            "median_umi_per_cell": [float(np.median(np.asarray(counts.sum(axis=1)).ravel()))],
            "counts_file": [str(count_file)],
        }
    )
    summary.to_csv(out_dir / f"{args.dataset_name}_glmpca_sample_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote sampled HVG counts to {count_file}")


if __name__ == "__main__":
    main()
