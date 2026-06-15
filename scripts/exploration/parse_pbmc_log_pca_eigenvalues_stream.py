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


def read_metadata(zip_file):
    proc = open_zip_member(zip_file, "cell_metadata_1M_PBMC.csv")
    metadata = pd.read_csv(proc.stdout, usecols=["tscp_count"])
    status = proc.wait()
    if status != 0:
        raise RuntimeError(f"unzip metadata stream exited with status {status}")
    return metadata


def finalize_range(rows, cols, data, range_start, range_end, gene_count, lib_size, accum):
    n_cells = range_end - range_start + 1
    if len(data) == 0:
        x = sp.csr_matrix((n_cells, gene_count), dtype=np.float64)
    else:
        x = sp.coo_matrix(
            (data, (rows, cols)),
            shape=(n_cells, gene_count),
            dtype=np.float64,
        ).tocsr()
        x.sum_duplicates()

    libs = lib_size[(range_start - 1):range_end].astype(np.float64, copy=True)
    libs[libs <= 0] = 1.0
    x = sp.diags(1e4 / libs) @ x
    x.data = np.log1p(x.data)

    accum["sum"] += np.asarray(x.sum(axis=0)).ravel()
    accum["xtx"] += (x.T @ x).toarray()
    accum["n"] += n_cells


def covariance_eigenvalues(accum):
    n = accum["n"]
    mean = accum["sum"] / n
    cov = (accum["xtx"] - n * np.outer(mean, mean)) / (n - 1)
    cov = (cov + cov.T) / 2
    eigvals = np.linalg.eigvalsh(cov)[::-1]
    eigvals[(eigvals < 0) & (eigvals > -1e-8)] = 0
    return eigvals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_file")
    parser.add_argument("selected_features_csv")
    parser.add_argument("out_dir")
    parser.add_argument("--dataset-name", default="parse_pbmc_log_pca")
    parser.add_argument("--cell-counts", default="5000,10000,25000,50000,100000,250000,500000,all")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(args.selected_features_csv)
    if "matrix_col" not in selected.columns:
        raise ValueError("selected_features_csv must contain a 1-based matrix_col column")
    selected_gene_ids = selected["matrix_col"].astype(int).to_numpy()
    selected_map = {int(gene): i for i, gene in enumerate(selected_gene_ids)}
    gene_count = len(selected_gene_ids)

    metadata = read_metadata(args.zip_file)
    lib_size = metadata["tscp_count"].to_numpy(dtype=np.float64)

    proc = open_zip_member(args.zip_file, "DGE_1M_PBMC.mtx")
    header = proc.stdout.readline().strip()
    if not header.startswith("%%MatrixMarket"):
        raise ValueError(f"Unexpected MatrixMarket header: {header}")
    line = proc.stdout.readline().strip()
    while line.startswith("%"):
        line = proc.stdout.readline().strip()
    n_total_cells, n_total_genes, n_total_nnz = map(int, line.split())

    checkpoints = []
    for item in args.cell_counts.split(","):
        item = item.strip().lower()
        if item in {"all", "inf"}:
            checkpoints.append(n_total_cells)
        elif item:
            checkpoints.append(min(int(item), n_total_cells))
    checkpoints = sorted(set(checkpoints))

    accum = {
        "sum": np.zeros(gene_count, dtype=np.float64),
        "xtx": np.zeros((gene_count, gene_count), dtype=np.float64),
        "n": 0,
    }

    eig_rows = []
    summary_rows = []
    checkpoint_idx = 0
    range_start = 1
    range_end = checkpoints[checkpoint_idx]
    rows = []
    cols = []
    data = []
    selected_nnz = 0

    def close_ranges_until(cell):
        nonlocal range_start, range_end, rows, cols, data, checkpoint_idx
        while checkpoint_idx < len(checkpoints) and cell > range_end:
            finalize_range(rows, cols, data, range_start, range_end, gene_count, lib_size, accum)
            eigvals = covariance_eigenvalues(accum)
            n_cells = accum["n"]
            eig_rows.append(
                pd.DataFrame(
                    {
                        "dataset": args.dataset_name,
                        "n_cells": n_cells,
                        "n_genes": gene_count,
                        "rank": np.arange(1, len(eigvals) + 1),
                        "eigenvalue": eigvals,
                        "log10_eigenvalue": np.log10(np.maximum(eigvals, np.finfo(float).tiny)),
                    }
                )
            )
            summary_rows.append(
                {
                    "dataset": args.dataset_name,
                    "n_total_cells": n_total_cells,
                    "n_cells": n_cells,
                    "n_total_genes": n_total_genes,
                    "n_total_nnz": n_total_nnz,
                    "n_genes_used": gene_count,
                    "mean_umi_per_cell": float(np.mean(lib_size[:n_cells])),
                    "median_umi_per_cell": float(np.median(lib_size[:n_cells])),
                    "max_eigenvalue": float(eigvals[0]),
                    "median_eigenvalue": float(np.median(eigvals)),
                    "mean_eigenvalue": float(np.mean(eigvals)),
                }
            )
            print(f"Checkpoint n={n_cells}: lambda1={eigvals[0]:.4f}, median={np.median(eigvals):.4f}", flush=True)
            checkpoint_idx += 1
            if checkpoint_idx >= len(checkpoints):
                return
            range_start = range_end + 1
            range_end = checkpoints[checkpoint_idx]
            rows = []
            cols = []
            data = []

    for raw_line in proc.stdout:
        cell_s, gene_s, value_s = raw_line.split()
        cell = int(cell_s)
        if checkpoint_idx >= len(checkpoints):
            break
        close_ranges_until(cell)
        if checkpoint_idx >= len(checkpoints):
            break
        gene = int(gene_s)
        mapped_gene = selected_map.get(gene)
        if mapped_gene is None:
            continue
        rows.append(cell - range_start)
        cols.append(mapped_gene)
        data.append(float(value_s))
        selected_nnz += 1

    if checkpoint_idx >= len(checkpoints):
        proc.terminate()
        proc.wait()
    else:
        status = proc.wait()
        if status != 0:
            raise RuntimeError(f"unzip DGE stream exited with status {status}")

    close_ranges_until(n_total_cells + 1)

    eig_df = pd.concat(eig_rows, ignore_index=True)
    summary_df = pd.DataFrame(summary_rows)

    eig_file = out_dir / f"{args.dataset_name}_log_pca_eigenvalues_by_n.csv"
    summary_file = out_dir / f"{args.dataset_name}_log_pca_eigenvalues_summary_by_n.csv"
    eig_df.to_csv(eig_file, index=False)
    summary_df.to_csv(summary_file, index=False)

    selected.to_csv(out_dir / f"{args.dataset_name}_selected_features.csv", index=False)
    print(summary_df.to_string(index=False))
    print(f"Selected-gene nnz streamed: {selected_nnz}")
    print(f"Wrote eigenvalues to {eig_file}")


if __name__ == "__main__":
    main()
