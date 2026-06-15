#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reusable Scanpy preprocessing + PCA template for single-cell count matrices."
    )
    parser.add_argument("--dataset-name", required=True, help="Short dataset name used in output files.")
    parser.add_argument(
        "--input-type",
        required=True,
        choices=["h5ad", "10x_mtx", "10x_h5", "csv", "tsv"],
        help="Input format.",
    )
    parser.add_argument("--input-path", required=True, help="Input file or directory path.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--batch-key", default=None, help="Optional obs column used for coloring plots.")
    parser.add_argument("--celltype-key", default=None, help="Optional obs column used for coloring plots.")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--no-umap", action="store_true", help="Skip neighbors/UMAP.")
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Optional random cell cap for quick template runs on very large datasets.",
    )
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def read_input(input_type: str, input_path: Path) -> sc.AnnData:
    if input_type == "h5ad":
        return sc.read_h5ad(input_path)
    if input_type == "10x_mtx":
        return sc.read_10x_mtx(input_path, var_names="gene_symbols", cache=False)
    if input_type == "10x_h5":
        return sc.read_10x_h5(input_path)
    if input_type in {"csv", "tsv"}:
        sep = "," if input_type == "csv" else "\t"
        df = pd.read_csv(input_path, sep=sep, index_col=0)
        return sc.AnnData(df)
    raise ValueError(f"Unsupported input type: {input_type}")


def first_present(adata: sc.AnnData, keys: list[str | None]) -> list[str]:
    return [key for key in keys if key and key in adata.obs.columns]


def save_plot(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)
    sc.settings.verbosity = 2
    sc.settings.figdir = str(out_dir)

    adata = read_input(args.input_type, input_path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    if args.max_cells is not None and adata.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)
        keep = np.sort(rng.choice(adata.n_obs, args.max_cells, replace=False))
        adata = adata[keep].copy()

    adata.layers["counts"] = adata.X.copy()
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    initial = {"cells": int(adata.n_obs), "features": int(adata.n_vars)}
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    filtered = {"cells": int(adata.n_obs), "features": int(adata.n_vars)}

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    n_hvg = min(args.n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata_hvg, max_value=10)
    n_pcs = min(args.n_pcs, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
    if n_pcs < 2:
        raise ValueError(f"Not enough cells/features after filtering for PCA: n_pcs={n_pcs}")
    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack", random_state=args.seed)

    if not args.no_umap:
        sc.pp.neighbors(adata_hvg, n_neighbors=args.neighbors, n_pcs=n_pcs, random_state=args.seed)
        sc.tl.umap(adata_hvg, random_state=args.seed)

    colors = first_present(adata_hvg, [args.celltype_key, args.batch_key])

    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=n_pcs, log=True, show=False)
    save_plot(out_dir / "pca_variance_ratio.png")

    if colors:
        sc.pl.pca(adata_hvg, color=colors, show=False)
    else:
        sc.pl.pca(adata_hvg, show=False)
    save_plot(out_dir / "pca_scatter.png")

    if "X_umap" in adata_hvg.obsm:
        if colors:
            sc.pl.umap(adata_hvg, color=colors, show=False)
        else:
            sc.pl.umap(adata_hvg, show=False)
        save_plot(out_dir / "umap.png")

    pca = adata_hvg.uns["pca"]
    variance_df = pd.DataFrame(
        {
            "pc": np.arange(1, len(pca["variance"]) + 1),
            "variance": pca["variance"],
            "variance_ratio": pca["variance_ratio"],
        }
    )
    variance_df.to_csv(out_dir / "pca_variance.csv", index=False)

    obs_summary = adata_hvg.obs.describe(include="all").transpose()
    obs_summary.to_csv(out_dir / "obs_summary.csv")

    metadata = {
        "dataset_name": args.dataset_name,
        "input_type": args.input_type,
        "input_path": str(input_path),
        "initial": initial,
        "filtered": filtered,
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "target_sum": args.target_sum,
        "n_hvg": int(adata_hvg.n_vars),
        "n_pcs": int(n_pcs),
        "neighbors": None if args.no_umap else args.neighbors,
        "batch_key": args.batch_key,
        "celltype_key": args.celltype_key,
        "seed": args.seed,
    }
    (out_dir / "processing_metadata.json").write_text(json.dumps(metadata, indent=2))

    adata_hvg.write_h5ad(out_dir / f"{args.dataset_name}_processed_pca.h5ad", compression="gzip")
    print(json.dumps(metadata, indent=2))
    print(f"Wrote reusable processed AnnData and PCA diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
