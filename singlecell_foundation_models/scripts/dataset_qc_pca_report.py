#!/usr/bin/env python
from __future__ import annotations

import argparse
import base64
import gzip
import html
import io
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a reusable HTML QC/PCA/MP/UMAP report for one single-cell dataset."
    )
    parser.add_argument("--dataset-name", required=True)
    parser.add_argument("--input-type", required=True, choices=["h5ad", "10x_mtx", "10x_h5", "csv", "tsv", "parse_zip"])
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--technology", default="unknown")
    parser.add_argument("--assay", default="scRNA-seq")
    parser.add_argument("--celltype-key", default=None)
    parser.add_argument("--batch-key", default=None)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument(
        "--embedding-pcs",
        type=int,
        default=30,
        help="Number of leading PCs used for neighbors, clustering, and UMAP.",
    )
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--transpose-table", action="store_true", help="Transpose csv/tsv input when rows are genes and columns are cells.")
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def choose_indices(n: int, max_cells: int | None, seed: int) -> np.ndarray:
    if max_cells is None or n <= max_cells:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_cells, replace=False))


def read_10x_mtx_sample(input_path: Path, max_cells: int | None, seed: int) -> sc.AnnData:
    barcodes = read_lines(input_path / "barcodes.tsv.gz")
    features_df = pd.read_csv(input_path / "features.tsv.gz", sep="\t", header=None)
    gene_names = features_df.iloc[:, 1].astype(str).to_numpy()
    selected = choose_indices(len(barcodes), max_cells, seed)
    selected_map = {old_idx + 1: new_idx for new_idx, old_idx in enumerate(selected)}

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    with gzip.open(input_path / "matrix.mtx.gz", "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_genes, n_cells, _ = map(int, line.split())
            break
        for line in handle:
            gene_idx, cell_idx, value = line.split()
            mapped = selected_map.get(int(cell_idx))
            if mapped is None:
                continue
            rows.append(mapped)
            cols.append(int(gene_idx) - 1)
            data.append(float(value))
    matrix = sparse.csr_matrix((data, (rows, cols)), shape=(len(selected), n_genes))
    adata = sc.AnnData(matrix)
    adata.obs_names = [barcodes[i] for i in selected]
    adata.var_names = gene_names
    return adata


def read_parse_zip_sample(input_path: Path, max_cells: int | None, seed: int) -> sc.AnnData:
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
            n_cells, n_genes, _ = map(int, line.split())
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
    adata = sc.AnnData(matrix)
    adata.obs = cell_meta.iloc[selected].copy()
    adata.obs_names = adata.obs["bc_wells"].astype(str).to_numpy()
    adata.var_names = genes["gene_name"].astype(str).to_numpy()
    return adata


def read_input(input_type: str, input_path: Path, max_cells: int | None, seed: int, transpose_table: bool) -> sc.AnnData:
    if input_type == "h5ad":
        if max_cells is not None:
            backed = sc.read_h5ad(input_path, backed="r")
            selected = choose_indices(backed.n_obs, max_cells, seed)
            adata = backed[selected, :].to_memory()
            backed.file.close()
            return adata
        return sc.read_h5ad(input_path)
    if input_type == "10x_mtx":
        if max_cells is not None:
            return read_10x_mtx_sample(input_path, max_cells, seed)
        return sc.read_10x_mtx(input_path, var_names="gene_symbols", cache=False)
    if input_type == "10x_h5":
        adata = sc.read_10x_h5(input_path)
        if "feature_types" in adata.var.columns:
            gene_mask = adata.var["feature_types"].astype(str).eq("Gene Expression").to_numpy()
            if gene_mask.any():
                adata = adata[:, gene_mask].copy()
        return adata
    if input_type in {"csv", "tsv"}:
        sep = "," if input_type == "csv" else "\t"
        df = pd.read_csv(input_path, sep=sep, index_col=0)
        if transpose_table:
            df = df.transpose()
        return sc.AnnData(df)
    if input_type == "parse_zip":
        return read_parse_zip_sample(input_path, max_cells, seed)
    raise ValueError(f"Unsupported input type: {input_type}")


def matrix_sum_per_row(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.sum(axis=1)).ravel()
    return np.asarray(x).sum(axis=1)


def matrix_nnz_per_row(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.diff(x.tocsr().indptr)
    return (np.asarray(x) > 0).sum(axis=1)


def matrix_min(x) -> float:
    if sparse.issparse(x):
        data_min = float(x.data.min()) if x.nnz else 0.0
        return min(0.0, data_min)
    return float(np.nanmin(np.asarray(x)))


def select_top_variable_genes(adata: sc.AnnData, n_hvg: int) -> None:
    x = adata.X
    if sparse.issparse(x):
        means = np.asarray(x.mean(axis=0)).ravel()
        means_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
        variances = means_sq - np.square(means)
    else:
        variances = np.nanvar(np.asarray(x), axis=0)
    variances = np.nan_to_num(variances, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
    n_pick = min(n_hvg, adata.n_vars)
    top = np.argsort(variances)[::-1][:n_pick]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top] = True
    adata.var["hvg_variance"] = variances
    adata.var["highly_variable"] = mask


def fmt_num(value: float | int | None, digits: int = 2) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    if abs(float(value)) >= 1e9:
        return f"{float(value) / 1e9:.{digits}f}B"
    if abs(float(value)) >= 1e6:
        return f"{float(value) / 1e6:.{digits}f}M"
    if abs(float(value)) >= 1e3:
        return f"{float(value):,.0f}"
    if float(value).is_integer():
        return f"{int(value)}"
    return f"{float(value):,.{digits}f}"


def mp_density(x: np.ndarray, gamma: float, sigma2: float) -> np.ndarray:
    lower = sigma2 * (1 - np.sqrt(gamma)) ** 2
    upper = sigma2 * (1 + np.sqrt(gamma)) ** 2
    out = np.zeros_like(x, dtype=float)
    inside = (x >= lower) & (x <= upper) & (x > 0)
    out[inside] = np.sqrt((upper - x[inside]) * (x[inside] - lower)) / (
        2 * np.pi * gamma * sigma2 * x[inside]
    )
    return out


def find_scree_knee(eigvals: np.ndarray, max_rank: int = 100) -> int:
    vals = np.sort(eigvals[np.isfinite(eigvals)])[::-1]
    k = min(max_rank, len(vals) - 1)
    if k < 3:
        return 0
    y = np.log10(np.maximum(vals[:k], np.finfo(float).tiny))
    x = np.arange(1, k + 1, dtype=float)
    x = (x - x.min()) / (x.max() - x.min())
    y = (y - y.min()) / (y.max() - y.min())
    x1, y1 = x[0], y[0]
    x2, y2 = x[-1], y[-1]
    dist = np.abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / np.sqrt(
        (y2 - y1) ** 2 + (x2 - x1) ** 2
    )
    return int(np.argmax(dist) + 1)


def estimate_mp_after_spikes(eigvals: np.ndarray, n_cells: int, k_spikes: int) -> dict[str, float | int]:
    vals = np.sort(eigvals[np.isfinite(eigvals)])[::-1]
    noise = vals[k_spikes:]
    noise = noise[noise > 1e-10]
    p_eff = len(noise)
    gamma = p_eff / n_cells
    sigma2 = float(np.mean(noise)) if p_eff else float("nan")
    lower = sigma2 * (1 - np.sqrt(gamma)) ** 2
    upper = sigma2 * (1 + np.sqrt(gamma)) ** 2
    return {
        "n_cells": int(n_cells),
        "n_eigenvalues": int(len(vals)),
        "k_spikes": int(k_spikes),
        "n_noise_positive": int(p_eff),
        "gamma": float(gamma),
        "sigma2": float(sigma2),
        "mp_lower": float(lower),
        "mp_upper": float(upper),
        "n_outliers_above_mp": int(np.sum(vals > upper)),
        "lambda1_over_mp_upper": float(vals[0] / upper) if upper > 0 and len(vals) else float("nan"),
    }


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=190, bbox_inches="tight")
    plt.close()


def img_tag(path: Path, alt: str) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" alt="{html.escape(alt)}">'


def table_html(rows: list[tuple[str, str]]) -> str:
    body = "\n".join(
        f"<tr><th>{html.escape(k)}</th><td>{v}</td></tr>" for k, v in rows
    )
    return f"<table>{body}</table>"


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    np.random.seed(args.seed)

    adata = read_input(args.input_type, input_path, args.max_cells, args.seed, args.transpose_table)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    if args.max_cells is not None and adata.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)
        keep = np.sort(rng.choice(adata.n_obs, args.max_cells, replace=False))
        adata = adata[keep].copy()

    raw_counts = adata.X.copy()
    umi = matrix_sum_per_row(raw_counts).astype(float)
    detected = matrix_nnz_per_row(raw_counts).astype(float)
    count_like = matrix_min(raw_counts) >= 0
    adata.obs["total_counts_report"] = umi
    adata.obs["n_features_report"] = detected

    qc_rows = [
        ("Dataset", html.escape(args.dataset_name)),
        ("Technology", html.escape(args.technology)),
        ("Assay", html.escape(args.assay)),
        ("Input", f"<code>{html.escape(str(input_path))}</code>"),
        ("Cells before filtering", fmt_num(adata.n_obs)),
        ("Features before filtering", fmt_num(adata.n_vars)),
        ("Total UMIs/counts", fmt_num(float(np.sum(umi)))),
        ("Median UMI/cell", fmt_num(float(np.median(umi)))),
        ("Mean UMI/cell", fmt_num(float(np.mean(umi)))),
        ("Median detected features/cell", fmt_num(float(np.median(detected)))),
    ]

    qc_vector = umi if count_like else np.abs(umi)
    qc_label = "UMI/counts per cell" if count_like else "absolute expression sum per cell"
    rank = np.arange(1, len(qc_vector) + 1)
    sorted_umi = np.sort(np.maximum(qc_vector, 1))[::-1]
    plt.figure(figsize=(7, 4.8))
    plt.plot(rank, sorted_umi, color="#2454a6", linewidth=1.0)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Cell barcode rank")
    plt.ylabel(qc_label)
    plt.title(f"{args.dataset_name}: UMI rank plot")
    plt.grid(alpha=0.25)
    rank_plot = out_dir / "umi_rank_plot.png"
    savefig(rank_plot)

    plt.figure(figsize=(7, 4.8))
    plt.hist(np.log10(np.maximum(qc_vector, 1)), bins=80, color="#18745a", alpha=0.85)
    plt.xlabel(f"log10({qc_label})")
    plt.ylabel("Cells")
    plt.title(f"{args.dataset_name}: UMI distribution")
    umi_hist = out_dir / "umi_distribution.png"
    savefig(umi_hist)

    if count_like:
        sc.pp.filter_cells(adata, min_genes=args.min_genes)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
    filtered_cells, filtered_features = adata.n_obs, adata.n_vars

    adata.layers["counts"] = adata.X.copy()
    if count_like:
        sc.pp.normalize_total(adata, target_sum=args.target_sum)
        sc.pp.log1p(adata)

    n_hvg = min(args.n_hvg, adata.n_vars)
    select_top_variable_genes(adata, n_hvg)
    hvg_flavor = "top_variance_log_normalized" if count_like else "top_variance_input_matrix"
    hvg_fallback_reason = None
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata_hvg, max_value=10)

    full_pcs = min(adata_hvg.n_vars, adata_hvg.n_obs - 1)
    embedding_pcs = min(args.embedding_pcs, full_pcs)
    sc.tl.pca(adata_hvg, n_comps=full_pcs, svd_solver="randomized", random_state=args.seed)
    eigvals = np.asarray(adata_hvg.uns["pca"]["variance"], dtype=float)
    eig_df = pd.DataFrame(
        {
            "rank": np.arange(1, len(eigvals) + 1),
            "eigenvalue": eigvals,
            "log10_eigenvalue": np.log10(np.maximum(eigvals, np.finfo(float).tiny)),
        }
    )
    eig_df.to_csv(out_dir / "pca_eigenvalues.csv", index=False)

    k_spikes = find_scree_knee(eigvals)
    mp_fit = estimate_mp_after_spikes(eigvals, n_cells=adata_hvg.n_obs, k_spikes=k_spikes)
    pd.DataFrame([mp_fit]).to_csv(out_dir / "mp_fit_parameters.csv", index=False)

    x = np.linspace(max(mp_fit["mp_lower"], 1e-12), mp_fit["mp_upper"], 1000)
    curve = pd.DataFrame({"eigenvalue": x, "density": mp_density(x, mp_fit["gamma"], mp_fit["sigma2"])})
    curve.to_csv(out_dir / "mp_fit_curve.csv", index=False)

    plt.figure(figsize=(7, 4.8))
    plt.plot(eig_df["rank"], eig_df["eigenvalue"], marker="o", markersize=3, color="#17202a")
    plt.axvline(k_spikes, color="#9f2a3b", linestyle="--", label=f"scree knee = {k_spikes}")
    plt.axhline(mp_fit["mp_upper"], color="#2454a6", linestyle="--", label="MP upper edge")
    plt.yscale("log")
    plt.xlabel("PC rank")
    plt.ylabel("Eigenvalue")
    plt.title(f"{args.dataset_name}: PCA eigenvalues")
    plt.legend(frameon=False)
    scree_plot = out_dir / "pca_eigenvalue_scree.png"
    savefig(scree_plot)

    bulk = eigvals[(eigvals > 1e-10) & (eigvals <= max(mp_fit["mp_upper"] * 1.25, mp_fit["mp_upper"]))]
    plt.figure(figsize=(7, 4.8))
    plt.hist(bulk, bins=40, density=True, color="#d9e1e8", edgecolor="white")
    plt.plot(curve["eigenvalue"], curve["density"], color="#2454a6", linewidth=2)
    plt.axvline(mp_fit["mp_upper"], color="#9f2a3b", linestyle="--", label="MP upper edge")
    plt.xlabel("PCA eigenvalue")
    plt.ylabel("Density")
    plt.title(f"{args.dataset_name}: MP fit to residual eigenvalue bulk")
    plt.legend(frameon=False)
    mp_plot = out_dir / "pca_mp_bulk_fit.png"
    savefig(mp_plot)

    sc.pp.neighbors(adata_hvg, n_neighbors=args.neighbors, n_pcs=embedding_pcs, random_state=args.seed)
    sc.tl.leiden(
        adata_hvg,
        resolution=args.leiden_resolution,
        key_added="leiden",
        random_state=args.seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    sc.tl.umap(adata_hvg, random_state=args.seed)
    color_keys = ["leiden"]
    color_keys.extend([k for k in [args.celltype_key, args.batch_key] if k and k in adata_hvg.obs.columns])
    sc.pl.umap(adata_hvg, color=color_keys[:2] if color_keys else None, show=False)
    umap_plot = out_dir / "umap.png"
    savefig(umap_plot)

    adata.obs["leiden"] = adata_hvg.obs["leiden"].reindex(adata.obs_names).astype(str)

    adata_hvg.write_h5ad(out_dir / f"{args.dataset_name}_processed_pca.h5ad", compression="gzip")

    metadata = {
        "dataset_name": args.dataset_name,
        "technology": args.technology,
        "assay": args.assay,
        "input_type": args.input_type,
        "input_path": str(input_path),
        "cells_before_filtering": int(len(umi)),
        "features_before_filtering": int(raw_counts.shape[1]),
        "cells_after_filtering": int(filtered_cells),
        "features_after_filtering": int(filtered_features),
        "total_umis": float(np.sum(umi)),
        "median_umi_per_cell": float(np.median(umi)),
        "mean_umi_per_cell": float(np.mean(umi)),
        "median_features_per_cell": float(np.median(detected)),
        "count_like_matrix": bool(count_like),
        "n_hvg": int(adata_hvg.n_vars),
        "hvg_flavor": hvg_flavor,
        "hvg_fallback_reason": hvg_fallback_reason,
        "full_pca_components": int(full_pcs),
        "embedding_pcs": int(embedding_pcs),
        "mp_fit": mp_fit,
        "leiden_resolution": args.leiden_resolution,
        "n_leiden_clusters": int(adata_hvg.obs["leiden"].nunique()),
    }
    (out_dir / "report_metrics.json").write_text(json.dumps(metadata, indent=2))

    qc_rows.extend(
        [
            ("Cells after filtering", fmt_num(filtered_cells)),
            ("Features after filtering", fmt_num(filtered_features)),
            ("Count-like matrix", "yes" if count_like else "no; treated as transformed expression"),
            ("Highly variable genes", fmt_num(adata_hvg.n_vars)),
            ("HVG selection", html.escape(hvg_flavor)),
            ("Full PCA components stored", fmt_num(full_pcs)),
            ("PCs used for neighbors/UMAP", fmt_num(embedding_pcs)),
            ("Leiden clusters", fmt_num(adata_hvg.obs["leiden"].nunique())),
            ("MP spike cutoff rank", fmt_num(k_spikes)),
            ("MP upper edge", fmt_num(mp_fit["mp_upper"], digits=4)),
            ("Eigenvalue outliers above MP edge", fmt_num(mp_fit["n_outliers_above_mp"])),
            ("Top eigenvalue / MP edge", fmt_num(mp_fit["lambda1_over_mp_upper"], digits=2)),
        ]
    )

    html_out = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(args.dataset_name)} QC + PCA Report</title>
  <style>
    body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17202a; background: #f7f9fb; }}
    header {{ padding: 28px 34px; background: #fff; border-bottom: 1px solid #d9e1e8; }}
    main {{ max-width: 1220px; margin: 0 auto; padding: 24px; }}
    h1 {{ margin: 0 0 6px; font-size: 34px; }}
    h2 {{ margin: 0 0 12px; font-size: 22px; }}
    p {{ color: #667085; margin: 0; }}
    code {{ background: #eef3f7; border: 1px solid #dbe5ec; padding: 1px 5px; border-radius: 5px; }}
    .grid {{ display: grid; grid-template-columns: repeat(12, 1fr); gap: 16px; }}
    .card {{ grid-column: span 6; background: #fff; border: 1px solid #d9e1e8; border-radius: 8px; padding: 16px; }}
    .full {{ grid-column: 1 / -1; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ text-align: left; padding: 9px 10px; border-bottom: 1px solid #d9e1e8; vertical-align: top; }}
    th {{ color: #344054; background: #eef3f7; width: 260px; }}
    img {{ width: 100%; height: auto; border: 1px solid #d9e1e8; border-radius: 8px; background: white; }}
    .metric-row {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin: 16px 0; }}
    .metric {{ background: #fff; border: 1px solid #d9e1e8; border-radius: 8px; padding: 14px; }}
    .metric strong {{ display: block; font-size: 24px; }}
    .metric span {{ color: #667085; font-size: 13px; }}
    @media (max-width: 900px) {{ .card {{ grid-column: 1 / -1; }} .metric-row {{ grid-template-columns: 1fr 1fr; }} }}
  </style>
</head>
<body>
  <header>
    <h1>{html.escape(args.dataset_name)} QC + PCA Report</h1>
    <p>{html.escape(args.technology)} · {html.escape(args.assay)} · generated from <code>{html.escape(str(input_path))}</code></p>
  </header>
  <main>
    <div class="metric-row">
      <div class="metric"><strong>{fmt_num(len(umi))}</strong><span>cells before filtering</span></div>
      <div class="metric"><strong>{fmt_num(raw_counts.shape[1])}</strong><span>features before filtering</span></div>
      <div class="metric"><strong>{fmt_num(float(np.sum(umi)))}</strong><span>total UMIs/counts</span></div>
      <div class="metric"><strong>{fmt_num(mp_fit["n_outliers_above_mp"])}</strong><span>PC eigenvalue outliers</span></div>
    </div>
    <section class="grid">
      <article class="card full">
        <h2>Basic QC</h2>
        {table_html(qc_rows)}
      </article>
      <article class="card">
        <h2>UMI Rank Plot</h2>
        {img_tag(rank_plot, "UMI rank plot")}
      </article>
      <article class="card">
        <h2>UMI Distribution</h2>
        {img_tag(umi_hist, "UMI distribution")}
      </article>
      <article class="card">
        <h2>PCA Eigenvalues</h2>
        {img_tag(scree_plot, "PCA eigenvalue scree plot")}
      </article>
      <article class="card">
        <h2>MP Bulk Fit</h2>
        {img_tag(mp_plot, "MP fit to PCA eigenvalue bulk")}
      </article>
      <article class="card full">
        <h2>UMAP</h2>
        {img_tag(umap_plot, "UMAP")}
      </article>
    </section>
  </main>
</body>
</html>
"""
    (out_dir / "index.html").write_text(html_out)
    print(f"Wrote report: {out_dir / 'index.html'}")


if __name__ == "__main__":
    main()
