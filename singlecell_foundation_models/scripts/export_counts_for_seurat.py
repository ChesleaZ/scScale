#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import io
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io as scipy_io
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage count-like matrices as 10x-style mtx files for the R/Seurat report pipeline."
    )
    parser.add_argument("--input-type", required=True, choices=["h5ad", "parse_zip"])
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def choose_indices(n: int, max_cells: int | None, seed: int) -> np.ndarray:
    if max_cells is None or n <= max_cells:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_cells, replace=False))


def read_h5ad(input_path: Path, max_cells: int | None, seed: int) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    if max_cells is not None:
        backed = sc.read_h5ad(input_path, backed="r")
        selected = choose_indices(backed.n_obs, max_cells, seed)
        adata = backed[selected, :].to_memory()
        backed.file.close()
    else:
        adata = sc.read_h5ad(input_path)

    if adata.raw is not None:
        raw = adata.raw.to_adata()
        matrix = raw.X
        barcodes = raw.obs_names.astype(str).tolist()
        gene_names = raw.var_names.astype(str).tolist()
    else:
        if "feature_types" in adata.var.columns:
            gene_mask = adata.var["feature_types"].astype(str).eq("Gene Expression").to_numpy()
            if gene_mask.any():
                adata = adata[:, gene_mask].copy()
        matrix = adata.X
        barcodes = adata.obs_names.astype(str).tolist()
        gene_names = adata.var_names.astype(str).tolist()

    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    return matrix.tocsr(), barcodes, gene_names


def read_parse_zip(input_path: Path, max_cells: int | None, seed: int) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    cell_meta = pd.read_csv(
        io.BytesIO(subprocess.check_output(["unzip", "-p", str(input_path), "cell_metadata_1M_PBMC.csv"]))
    )
    genes = pd.read_csv(
        io.BytesIO(subprocess.check_output(["unzip", "-p", str(input_path), "all_genes_1M_PBMC.csv"]))
    )
    selected = choose_indices(len(cell_meta), max_cells, seed)
    selected_map = {old_idx + 1: new_idx for new_idx, old_idx in enumerate(selected)}

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    proc = subprocess.Popen(
        ["unzip", "-p", str(input_path), "DGE_1M_PBMC.mtx"],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            if line.startswith("%"):
                continue
            _n_cells, n_genes, _nnz = map(int, line.split())
            break
        for line in proc.stdout:
            cell_idx, gene_idx, value = line.split()
            mapped = selected_map.get(int(cell_idx))
            if mapped is None:
                continue
            rows.append(mapped)
            cols.append(int(gene_idx) - 1)
            data.append(float(value))
    finally:
        proc.stdout.close()
        return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(f"unzip failed while streaming Parse matrix with exit code {return_code}")

    matrix = sparse.csr_matrix((data, (rows, cols)), shape=(len(selected), n_genes))
    barcodes = cell_meta.iloc[selected]["bc_wells"].astype(str).tolist()
    gene_names = genes["gene_name"].astype(str).tolist()
    return matrix, barcodes, gene_names


def write_10x(matrix_cells_by_genes: sparse.csr_matrix, barcodes: list[str], genes: list[str], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_genes_by_cells = matrix_cells_by_genes.transpose().tocoo()
    with gzip.open(out_dir / "matrix.mtx.gz", "wb") as handle:
        scipy_io.mmwrite(handle, matrix_genes_by_cells)
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as handle:
        for barcode in barcodes:
            handle.write(f"{barcode}\n")
    with gzip.open(out_dir / "features.tsv.gz", "wt") as handle:
        for gene in genes:
            safe_gene = str(gene)
            handle.write(f"{safe_gene}\t{safe_gene}\tGene Expression\n")


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_path)
    out_dir = Path(args.out_dir)
    if args.input_type == "h5ad":
        matrix, barcodes, genes = read_h5ad(input_path, args.max_cells, args.seed)
    elif args.input_type == "parse_zip":
        matrix, barcodes, genes = read_parse_zip(input_path, args.max_cells, args.seed)
    else:
        raise ValueError(args.input_type)
    write_10x(matrix, barcodes, genes, out_dir)
    print(f"Wrote staged 10x matrix to {out_dir}")


if __name__ == "__main__":
    main()
