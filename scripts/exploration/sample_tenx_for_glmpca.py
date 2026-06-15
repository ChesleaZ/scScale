import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_matrix_shape(matrix_file):
    with gzip.open(matrix_file, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_genes, n_cells, n_nz = map(int, line.split())
            return n_genes, n_cells, n_nz
    raise ValueError(f"No MatrixMarket shape line found in {matrix_file}")


def stream_sampled_cells(matrix_file, n_genes, sampled_cells):
    sampled_set = set(sampled_cells.tolist())
    sampled_map = {cell: i for i, cell in enumerate(sampled_cells.tolist())}
    rows = []
    cols = []
    data = []

    with gzip.open(matrix_file, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_genes_seen, n_cells_seen, _ = map(int, line.split())
            if n_genes_seen != n_genes:
                raise ValueError("Unexpected MatrixMarket dimensions while streaming")
            break

        for line in handle:
            gene, cell, value = line.split()
            cell = int(cell)
            if cell not in sampled_set:
                continue
            rows.append(int(gene) - 1)
            cols.append(sampled_map[cell])
            data.append(int(float(value)))

    counts = sp.coo_matrix(
        (data, (rows, cols)),
        shape=(n_genes, len(sampled_cells)),
        dtype=np.int32,
    ).tocsr()
    counts.eliminate_zeros()
    return counts


def normalize_log1p_gene_by_cell(counts, scale_factor=1e4):
    lib_size = np.asarray(counts.sum(axis=0)).ravel()
    lib_size[lib_size <= 0] = 1
    x = counts.tocsc(copy=True) @ sp.diags(scale_factor / lib_size)
    x = x.tocsr()
    x.data = np.log1p(x.data)
    return x


def sparse_row_var(x):
    mean = np.asarray(x.mean(axis=1)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=1)).ravel()
    var = mean_sq - mean**2
    var[var < 0] = 0
    return var


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tenx_dir")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", default="jurkat")
    parser.add_argument("--sample-cells", type=int, default=1000)
    parser.add_argument("--n-variable-genes", type=int, default=400)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    tenx_dir = Path(args.tenx_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_file = tenx_dir / "matrix.mtx.gz"
    features_file = tenx_dir / "features.tsv.gz"
    barcodes_file = tenx_dir / "barcodes.tsv.gz"

    n_genes, n_cells, n_nz = read_matrix_shape(matrix_file)
    rng = np.random.default_rng(args.seed)
    sample_n = min(args.sample_cells, n_cells)
    sampled_cells = np.sort(rng.choice(n_cells, size=sample_n, replace=False) + 1)

    counts = stream_sampled_cells(matrix_file, n_genes, sampled_cells)
    x = normalize_log1p_gene_by_cell(counts)
    gene_var = sparse_row_var(x)
    n_hvg = min(args.n_variable_genes, n_genes)
    gene_idx = np.sort(np.argsort(gene_var)[-n_hvg:])

    counts_hvg = counts[gene_idx, :].toarray()

    features = pd.read_csv(features_file, sep="\t", header=None)
    with gzip.open(barcodes_file, "rt") as handle:
        barcodes = [line.strip() for line in handle]

    selected_features = features.iloc[gene_idx].copy()
    selected_features.columns = ["gene_id", "gene_name", "feature_type"][: selected_features.shape[1]]
    selected_features.insert(0, "matrix_row", gene_idx + 1)
    selected_features.to_csv(out_dir / f"{args.dataset_name}_glmpca_selected_features.csv", index=False)

    selected_barcodes = pd.DataFrame(
        {
            "matrix_col": sampled_cells,
            "barcode": [barcodes[i - 1] for i in sampled_cells],
            "umi_count": np.asarray(counts.sum(axis=0)).ravel(),
        }
    )
    selected_barcodes.to_csv(out_dir / f"{args.dataset_name}_glmpca_selected_cells.csv", index=False)

    gene_labels = selected_features["gene_name"].fillna(selected_features["gene_id"]).astype(str)
    gene_labels = [f"{label}_{idx}" for label, idx in zip(gene_labels, gene_idx + 1)]
    cell_labels = [f"cell_{i}" for i in sampled_cells]
    count_df = pd.DataFrame(counts_hvg, index=gene_labels, columns=cell_labels)
    count_file = out_dir / f"{args.dataset_name}_glmpca_counts_hvg.csv"
    count_df.to_csv(count_file)

    summary = pd.DataFrame(
        {
            "dataset": [args.dataset_name],
            "tenx_dir": [str(tenx_dir)],
            "n_total_genes": [n_genes],
            "n_total_cells": [n_cells],
            "n_total_nnz": [n_nz],
            "sampled_cells": [sample_n],
            "sampled_nnz": [counts.nnz],
            "n_genes_used": [n_hvg],
            "mean_umi_per_cell": [float(np.asarray(counts.sum(axis=0)).mean())],
            "median_umi_per_cell": [float(np.median(np.asarray(counts.sum(axis=0)).ravel()))],
            "counts_file": [str(count_file)],
        }
    )
    summary.to_csv(out_dir / f"{args.dataset_name}_glmpca_sample_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote sampled HVG counts to {count_file}")


if __name__ == "__main__":
    main()
