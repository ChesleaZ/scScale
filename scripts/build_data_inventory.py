#!/usr/bin/env python3
"""Build a lightweight HTML/CSV inventory for the shared single-cell data."""

from __future__ import annotations

import csv
import gzip
import html
import json
import math
import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Iterable

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "shared_data"


def human_bytes(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(num_bytes)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{value:.1f} TB"


def fmt_int(value: float | int | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{float(value):,.1f}"


def fmt_compact(value: float | int | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    value = float(value)
    units = ["", "K", "M", "B", "T"]
    unit_idx = 0
    while abs(value) >= 1000 and unit_idx < len(units) - 1:
        value /= 1000
        unit_idx += 1
    if unit_idx == 0:
        return f"{int(round(value)):,}"
    return f"{value:.1f}{units[unit_idx]}"


def line_count(path: Path, gzip_text: bool = False) -> int:
    opener = gzip.open if gzip_text else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _ in handle)


def matrix_market_dims(handle: BinaryIO) -> tuple[int, int, int]:
    first = handle.readline()
    if not first.startswith(b"%%MatrixMarket"):
        raise ValueError("not a Matrix Market file")
    for line in handle:
        if line.startswith(b"%"):
            continue
        rows, cols, nnz = line.split()[:3]
        return int(rows), int(cols), int(nnz)
    raise ValueError("missing Matrix Market dimensions")


def stream_matrix_market_by_cell(path: Path) -> tuple[int, int, int, np.ndarray, np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if not first.startswith("%%MatrixMarket"):
            raise ValueError("not a Matrix Market file")
        for line in handle:
            if line.startswith("%"):
                continue
            rows, cols, nnz = (int(v) for v in line.split()[:3])
            break
        else:
            raise ValueError("missing Matrix Market dimensions")

        umi_by_cell = np.zeros(cols, dtype=np.float64)
        genes_by_cell = np.zeros(cols, dtype=np.float64)
        for line in handle:
            if not line.strip():
                continue
            _, col, value = line.split()[:3]
            cell_idx = int(col) - 1
            umi_by_cell[cell_idx] += float(value)
            genes_by_cell[cell_idx] += 1

    return rows, cols, nnz, umi_by_cell, genes_by_cell


def scan_tenx_matrix(dataset: str, rel_dir: str) -> dict[str, object]:
    base = DATA_DIR / rel_dir
    matrix_path = base / "sample_filtered_feature_bc_matrix" / "matrix.mtx.gz"
    barcodes_path = base / "sample_filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    features_path = base / "sample_filtered_feature_bc_matrix" / "features.tsv.gz"

    print(f"Scanning {dataset} 10x matrix...", flush=True)
    barcode_count = line_count(barcodes_path, gzip_text=True)
    feature_count = line_count(features_path, gzip_text=True)

    rows, cols, nnz, umi_by_cell, genes_by_cell = stream_matrix_market_by_cell(matrix_path)

    return summarize_dataset(
        dataset=dataset,
        rel_path=rel_dir,
        source_format="10x Matrix Market",
        n_cells=barcode_count,
        n_genes=feature_count,
        matrix_rows=rows,
        matrix_cols=cols,
        nnz=nnz,
        umi_by_cell=umi_by_cell,
        genes_by_cell=genes_by_cell,
        size_bytes=dir_size(DATA_DIR / rel_dir),
        notes="Cells from barcodes.tsv.gz; genes from features.tsv.gz; UMIs streamed from matrix.mtx.gz.",
    )


def h5_dataset_shape(h5: h5py.File) -> tuple[int, int]:
    x = h5["X"]
    if isinstance(x, h5py.Dataset):
        return tuple(int(v) for v in x.shape)
    shape = x.attrs.get("shape")
    if shape is None:
        raise ValueError("H5AD X group lacks shape attribute")
    return int(shape[0]), int(shape[1])


def read_obs_vector(h5: h5py.File, column: str) -> np.ndarray | None:
    if "obs" not in h5 or column not in h5["obs"]:
        return None
    obj = h5["obs"][column]
    if not isinstance(obj, h5py.Dataset):
        return None
    return np.asarray(obj[()])


def scan_h5ad(dataset: str, rel_path: str, umi_column: str = "UMI_count", gene_column: str | None = None) -> dict[str, object]:
    path = DATA_DIR / rel_path
    with h5py.File(path, "r") as h5:
        n_cells, n_genes = h5_dataset_shape(h5)
        umi_by_cell = read_obs_vector(h5, umi_column)
        genes_by_cell = read_obs_vector(h5, gene_column) if gene_column else None
        x = h5["X"]
        storage = "dense" if isinstance(x, h5py.Dataset) else str(x.attrs.get("encoding-type", "sparse"))
        nnz = None
        if isinstance(x, h5py.Group) and "data" in x:
            nnz = int(x["data"].shape[0])

    notes = []
    if umi_by_cell is not None:
        notes.append(f"Per-cell counts from obs/{umi_column}.")
    else:
        notes.append(f"Matrix shape and nonzero count read from H5AD; obs/{umi_column} was not present.")
    if genes_by_cell is not None and gene_column:
        notes.append(f"Genes/cell from obs/{gene_column}.")

    return summarize_dataset(
        dataset=dataset,
        rel_path=rel_path,
        source_format=f"H5AD ({storage})",
        n_cells=n_cells,
        n_genes=n_genes,
        matrix_rows=n_cells,
        matrix_cols=n_genes,
        nnz=nnz,
        umi_by_cell=umi_by_cell,
        genes_by_cell=genes_by_cell,
        size_bytes=path.stat().st_size,
        notes=" ".join(notes),
    )


def scan_tenx_h5(dataset: str, rel_path: str, notes: str) -> dict[str, object]:
    path = DATA_DIR / rel_path
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        rows, cols = (int(v) for v in matrix["shape"][()])
        indptr = np.asarray(matrix["indptr"][()], dtype=np.int64)
        data = np.asarray(matrix["data"][()], dtype=np.float64)
        nnz = int(data.size)
        umi_by_cell = np.add.reduceat(data, indptr[:-1])
        genes_by_cell = np.diff(indptr)

    return summarize_dataset(
        dataset=dataset,
        rel_path=rel_path,
        source_format="10x feature-barcode HDF5",
        n_cells=cols,
        n_genes=rows,
        matrix_rows=rows,
        matrix_cols=cols,
        nnz=nnz,
        umi_by_cell=umi_by_cell,
        genes_by_cell=genes_by_cell,
        size_bytes=path.stat().st_size,
        notes=notes,
    )


def scan_wide_count_table(dataset: str, rel_path: str, notes: str) -> dict[str, object]:
    path = DATA_DIR / rel_path
    print(f"Scanning {dataset} wide count table...", flush=True)
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n_cells = len(header)
        umi_by_cell = np.zeros(n_cells, dtype=np.float64)
        genes_by_cell = np.zeros(n_cells, dtype=np.float64)
        n_genes = 0
        nnz = 0
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= 1:
                continue
            values = np.fromiter((float(v or 0) for v in parts[1:]), dtype=np.float64, count=n_cells)
            umi_by_cell += values
            detected = values > 0
            genes_by_cell += detected
            nnz += int(np.count_nonzero(detected))
            n_genes += 1

    return summarize_dataset(
        dataset=dataset,
        rel_path=rel_path,
        source_format="Tab-delimited gene x cell count table",
        n_cells=n_cells,
        n_genes=n_genes,
        matrix_rows=n_genes,
        matrix_cols=n_cells,
        nnz=nnz,
        umi_by_cell=umi_by_cell,
        genes_by_cell=genes_by_cell,
        size_bytes=path.stat().st_size,
        notes=notes,
    )


def unzip_member_lines(zip_path: Path, member: str) -> Iterable[str]:
    proc = subprocess.Popen(
        ["unzip", "-p", str(zip_path), member],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            yield line
    finally:
        proc.stdout.close()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        return_code = proc.wait()
        if return_code != 0:
            raise RuntimeError(f"unzip failed for {member}: {stderr.strip()}")


def zip_member_matrix_dims(zip_path: Path, member: str) -> tuple[int, int, int]:
    proc = subprocess.Popen(
        ["unzip", "-p", str(zip_path), member],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.stdout is not None
    try:
        dims = matrix_market_dims(proc.stdout)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
    return dims


def scan_parse_zip() -> dict[str, object]:
    rel_path = "Parse_1M_PBMC_T1D/1M_PBMC_T1D_Parse.zip"
    zip_path = DATA_DIR / rel_path
    umi_values: list[float] = []
    gene_values: list[float] = []

    reader = csv.DictReader(unzip_member_lines(zip_path, "cell_metadata_1M_PBMC.csv"))
    for row in reader:
        umi_values.append(float(row["tscp_count"]))
        gene_values.append(float(row["gene_count"]))

    gene_rows = max(0, sum(1 for _ in unzip_member_lines(zip_path, "all_genes_1M_PBMC.csv")) - 1)
    matrix_rows, matrix_cols, nnz = zip_member_matrix_dims(zip_path, "DGE_1M_PBMC.mtx")

    return summarize_dataset(
        dataset="Parse 1M PBMC T1D",
        rel_path=rel_path,
        source_format="ZIP: Matrix Market + cell metadata",
        n_cells=len(umi_values),
        n_genes=gene_rows,
        matrix_rows=matrix_rows,
        matrix_cols=matrix_cols,
        nnz=nnz,
        umi_by_cell=np.asarray(umi_values, dtype=np.float64),
        genes_by_cell=np.asarray(gene_values, dtype=np.float64),
        size_bytes=zip_path.stat().st_size,
        notes="Per-cell UMIs from cell_metadata_1M_PBMC.csv tscp_count; genes/cell from gene_count; matrix dimensions from DGE_1M_PBMC.mtx header.",
    )


def summarize_dataset(
    *,
    dataset: str,
    rel_path: str,
    source_format: str,
    n_cells: int,
    n_genes: int,
    matrix_rows: int,
    matrix_cols: int,
    nnz: int | None,
    umi_by_cell: np.ndarray | None,
    genes_by_cell: np.ndarray | None,
    size_bytes: int,
    notes: str,
) -> dict[str, object]:
    umi_by_cell = None if umi_by_cell is None else np.asarray(umi_by_cell, dtype=np.float64)
    genes_by_cell = None if genes_by_cell is None else np.asarray(genes_by_cell, dtype=np.float64)

    total_umis = float(np.nansum(umi_by_cell)) if umi_by_cell is not None else None
    mean_umi = float(np.nanmean(umi_by_cell)) if umi_by_cell is not None and umi_by_cell.size else None
    median_umi = float(np.nanmedian(umi_by_cell)) if umi_by_cell is not None and umi_by_cell.size else None
    mean_genes = float(np.nanmean(genes_by_cell)) if genes_by_cell is not None and genes_by_cell.size else None
    median_genes = float(np.nanmedian(genes_by_cell)) if genes_by_cell is not None and genes_by_cell.size else None

    return {
        "dataset": dataset,
        "path": rel_path,
        "format": source_format,
        "size_bytes": size_bytes,
        "size": human_bytes(size_bytes),
        "cells": int(n_cells),
        "genes": int(n_genes),
        "matrix_rows": int(matrix_rows),
        "matrix_cols": int(matrix_cols),
        "nonzero_entries": int(nnz) if nnz is not None else None,
        "total_umis": total_umis,
        "mean_umis_per_cell": mean_umi,
        "median_umis_per_cell": median_umi,
        "mean_genes_per_cell": mean_genes,
        "median_genes_per_cell": median_genes,
        "notes": notes,
    }


def dir_size(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            total += (Path(root) / name).stat().st_size
    return total


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    fields = [
        "dataset",
        "path",
        "format",
        "size_bytes",
        "size",
        "cells",
        "genes",
        "matrix_rows",
        "matrix_cols",
        "nonzero_entries",
        "total_umis",
        "mean_umis_per_cell",
        "median_umis_per_cell",
        "mean_genes_per_cell",
        "median_genes_per_cell",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_html(rows: list[dict[str, object]], path: Path) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    total_size = sum(int(row["size_bytes"]) for row in rows)
    total_cells = sum(int(row["cells"]) for row in rows)
    total_umis = sum(float(row["total_umis"] or 0) for row in rows)
    max_cells = max(float(row["cells"]) for row in rows)
    headers = [
        ("dataset", "Dataset"),
        ("size", "Size"),
        ("cells", "Cells"),
        ("genes", "Genes"),
        ("median_umis_per_cell", "Median UMI/cell"),
        ("mean_umis_per_cell", "Mean UMI/cell"),
        ("median_genes_per_cell", "Median genes/cell"),
        ("mean_genes_per_cell", "Mean genes/cell"),
        ("total_umis", "Total UMIs"),
        ("nonzero_entries", "Nonzero"),
        ("path", "Path"),
        ("notes", "Notes"),
    ]
    body_rows = []
    for row in rows:
        search_text = " ".join(str(row.get(key, "")) for key, _ in headers)
        cells_share = 100 * float(row["cells"]) / max_cells if max_cells else 0
        cells = []
        for key, label in headers:
            value = row.get(key)
            if key == "dataset":
                cells.append(
                    "<td class=\"dataset-cell\">"
                    f"<strong>{html.escape(str(value))}</strong>"
                    f"<span>{html.escape(str(row.get('format', '')))}</span>"
                    "</td>"
                )
            elif key == "path":
                href = "../data/" + html.escape(str(value), quote=True)
                display = html.escape(str(value))
                cells.append(f'<td class="path-cell"><a href="{href}">{display}</a></td>')
            elif key in {"cells", "genes", "total_umis", "nonzero_entries"}:
                if key == "cells":
                    cells.append(
                        "<td class=\"num metric-cell\">"
                        f"<span>{fmt_int(value)}</span>"
                        f"<i style=\"--bar:{cells_share:.2f}%\"></i>"
                        "</td>"
                    )
                else:
                    cells.append(f"<td class=\"num\">{fmt_int(value)}</td>")
            elif key in {
                "median_umis_per_cell",
                "mean_umis_per_cell",
                "median_genes_per_cell",
                "mean_genes_per_cell",
            }:
                cells.append(f"<td class=\"num\">{fmt_float(value)}</td>")
            else:
                cells.append(f"<td>{html.escape(str(value) if value is not None else '')}</td>")
        body_rows.append(
            f'<tr data-search="{html.escape(search_text.lower(), quote=True)}">'
            + "".join(cells)
            + "</tr>"
        )

    header_html = "".join(f"<th>{html.escape(label)}</th>" for _, label in headers)
    rows_html = "\n".join(body_rows)
    path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>RMT_sc Shared Data Inventory</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #172126;
      --muted: #5f6c73;
      --subtle: #eef3f0;
      --line: #d4ded9;
      --line-strong: #aebeb8;
      --head: #e7f0ec;
      --accent: #0b6f68;
      --accent-soft: #d6ece7;
      --bg: #f7f9f8;
      --warn: #816b2e;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    main {{
      max-width: 1440px;
      margin: 0 auto;
      padding: 26px 24px 40px;
    }}
    header {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(260px, 420px);
      gap: 20px;
      align-items: end;
      margin-bottom: 18px;
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 28px;
      font-weight: 700;
      letter-spacing: 0;
    }}
    .meta {{
      margin: 0;
      color: var(--muted);
    }}
    .search {{
      display: flex;
      align-items: center;
      gap: 9px;
      min-height: 40px;
      padding: 0 12px;
      border: 1px solid var(--line-strong);
      background: white;
    }}
    .search svg {{
      flex: 0 0 auto;
      color: var(--muted);
    }}
    .search input {{
      width: 100%;
      min-width: 0;
      border: 0;
      outline: 0;
      font: inherit;
      color: var(--ink);
      background: transparent;
    }}
    .stats {{
      display: grid;
      grid-template-columns: repeat(4, minmax(150px, 1fr));
      border: 1px solid var(--line);
      background: white;
      margin-bottom: 18px;
    }}
    .stat {{
      padding: 13px 16px;
      border-right: 1px solid var(--line);
    }}
    .stat:last-child {{
      border-right: 0;
    }}
    .stat .label {{
      display: block;
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
    }}
    .stat .value {{
      display: block;
      margin-top: 3px;
      font-size: 22px;
      font-weight: 700;
      font-variant-numeric: tabular-nums;
    }}
    .table-wrap {{
      overflow-x: auto;
      border: 1px solid var(--line);
      background: white;
      box-shadow: 0 1px 2px rgba(23, 33, 38, 0.04);
    }}
    table {{
      width: 100%;
      min-width: 1280px;
      border-collapse: collapse;
    }}
    th, td {{
      padding: 11px 12px;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
      text-align: left;
    }}
    th {{
      position: sticky;
      top: 0;
      background: var(--head);
      font-size: 12px;
      text-transform: uppercase;
      color: #344349;
      z-index: 1;
      white-space: nowrap;
    }}
    tbody tr:hover {{
      background: #fbfdfc;
    }}
    tbody tr.is-hidden {{
      display: none;
    }}
    .dataset-cell {{
      min-width: 190px;
    }}
    .dataset-cell strong,
    .dataset-cell span {{
      display: block;
    }}
    .dataset-cell span {{
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
    }}
    .path-cell {{
      max-width: 300px;
      overflow-wrap: anywhere;
    }}
    td.num {{
      text-align: right;
      white-space: nowrap;
      font-variant-numeric: tabular-nums;
    }}
    .metric-cell span {{
      display: block;
      margin-bottom: 5px;
    }}
    .metric-cell i {{
      display: block;
      width: 100%;
      height: 5px;
      background: var(--subtle);
      overflow: hidden;
    }}
    .metric-cell i::before {{
      content: "";
      display: block;
      width: var(--bar);
      height: 100%;
      background: var(--accent);
    }}
    a {{
      color: var(--accent);
      text-decoration: none;
    }}
    a:hover {{
      text-decoration: underline;
    }}
    code {{
      padding: 1px 4px;
      background: var(--subtle);
      border: 1px solid var(--line);
      font-size: 12px;
    }}
    .footnote {{
      margin-top: 14px;
      color: var(--muted);
      max-width: 980px;
    }}
    .empty {{
      display: none;
      padding: 18px;
      border: 1px solid var(--line);
      border-top: 0;
      background: white;
      color: var(--warn);
    }}
    @media (max-width: 820px) {{
      main {{
        padding: 20px 14px 30px;
      }}
      header {{
        grid-template-columns: 1fr;
        align-items: stretch;
      }}
      .stats {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
      .stat:nth-child(2) {{
        border-right: 0;
      }}
      .stat:nth-child(-n+2) {{
        border-bottom: 1px solid var(--line);
      }}
    }}
  </style>
</head>
<body>
  <main>
    <header>
      <div>
        <h1>RMT_sc Shared Data Inventory</h1>
        <p class="meta">Generated {html.escape(generated)} from files under <code>data/</code>.</p>
      </div>
      <label class="search" aria-label="Search datasets">
        <svg aria-hidden="true" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
          <circle cx="11" cy="11" r="8"></circle>
          <path d="m21 21-4.3-4.3"></path>
        </svg>
        <input id="search" type="search" placeholder="Search dataset, format, path, or notes">
      </label>
    </header>
    <section class="stats" aria-label="Inventory summary">
      <div class="stat"><span class="label">Datasets</span><span class="value">{len(rows)}</span></div>
      <div class="stat"><span class="label">Cells</span><span class="value">{fmt_compact(total_cells)}</span></div>
      <div class="stat"><span class="label">UMIs</span><span class="value">{fmt_compact(total_umis)}</span></div>
      <div class="stat"><span class="label">Data size</span><span class="value">{human_bytes(total_size)}</span></div>
    </section>
    <div class="table-wrap">
      <table id="inventory">
        <thead><tr>{header_html}</tr></thead>
        <tbody>
{rows_html}
        </tbody>
      </table>
    </div>
    <div id="empty" class="empty">No datasets match this search.</div>
    <p class="footnote">UMI and genes-per-cell summaries use per-cell metadata when present. For 10x Matrix Market inputs, UMIs and detected genes per cell are computed by streaming the sparse matrix.</p>
  </main>
  <script>
    const search = document.querySelector("#search");
    const rows = Array.from(document.querySelectorAll("#inventory tbody tr"));
    const empty = document.querySelector("#empty");

    search.addEventListener("input", () => {{
      const terms = search.value.trim().toLowerCase().split(/\\s+/).filter(Boolean);
      let visible = 0;
      rows.forEach((row) => {{
        const haystack = row.dataset.search || "";
        const match = terms.every((term) => haystack.includes(term));
        row.classList.toggle("is-hidden", !match);
        if (match) visible += 1;
      }});
      empty.style.display = visible === 0 ? "block" : "none";
    }});
  </script>
</body>
</html>
""",
        encoding="utf-8",
    )


def load_cached_rows(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    try:
        rows = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return {str(row.get("dataset")): row for row in rows if isinstance(row, dict)}


def cached_or_scan(cache: dict[str, dict[str, object]], dataset: str, scanner) -> dict[str, object]:
    cached = cache.get(dataset)
    if cached is not None:
        print(f"Using cached inventory row for {dataset}...", flush=True)
        return cached
    return scanner()


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    cache = load_cached_rows(OUT_DIR / "data_inventory.json")
    rows = [
        cached_or_scan(cache, "Jurkat", lambda: scan_tenx_matrix("Jurkat", "Jurkat")),
        cached_or_scan(cache, "K562", lambda: scan_tenx_matrix("K562", "K562")),
        cached_or_scan(
            cache,
            "Norman-Weissman 2019",
            lambda: scan_h5ad("Norman-Weissman 2019", "NormanWeissman2019/NormanWeissman2019_filtered.h5ad", gene_column="ngenes"),
        ),
        cached_or_scan(cache, "Parse 1M PBMC T1D", scan_parse_zip),
        cached_or_scan(cache, "Replogle RPE1", lambda: scan_h5ad("Replogle RPE1", "Replogle_RPE1/rpe1_raw_singlecell_01.h5ad")),
        scan_h5ad("FedscGen PBMC", "FedscGen/PBMC.h5ad", umi_column="n_counts", gene_column="n_genes"),
        scan_h5ad("FedscGen Human Pancreas", "FedscGen/HumanPancreas.h5ad", umi_column="n_counts", gene_column="n_genes"),
        scan_h5ad(
            "BMMC Multiome GEX",
            "Bone_Marrow_BMMCs/multiome_gex_processed_training.h5ad",
            umi_column="n_counts",
            gene_column="n_genes",
        ),
        scan_h5ad(
            "BMMC Multiome ATAC",
            "Bone_Marrow_BMMCs/multiome_atac_processed_training.h5ad",
            umi_column="nCount_peaks",
            gene_column=None,
        ),
        scan_tenx_h5(
            "10x PBMC Multiome 3k",
            "10x_PBMC_Multiome_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5",
            "UMIs and detected features computed from the 10x feature-barcode HDF5 sparse matrix.",
        ),
        scan_h5ad(
            "Human Lung Cell Atlas core",
            "HLCA/HLCA_core_CELLxGENE_2025-11-08.h5ad",
            umi_column="n_counts",
            gene_column=None,
        ),
        scan_wide_count_table(
            "Smart-seq3 PBMC UMI counts",
            "SmartSeq3_E-MTAB-8735/HCA.UMIcounts.PBMC.txt",
            "High-depth Smart-seq3 non-10x PBMC UMI count table from E-MTAB-8735; cells from columns and genes from rows.",
        ),
    ]
    write_csv(rows, OUT_DIR / "data_inventory.csv")
    (OUT_DIR / "data_inventory.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    write_html(rows, OUT_DIR / "index.html")
    print(f"Wrote {OUT_DIR / 'index.html'}")


if __name__ == "__main__":
    main()
