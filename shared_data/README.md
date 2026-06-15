# Shared Data Inventory

This directory is the human-facing index for the shared `RMT_sc` datasets.

- `index.html`: browser-friendly table of dataset metadata.
- `data_inventory.csv`: compact tabular copy of the same metadata.
- `data_inventory.json`: machine-readable copy of the same metadata.

Regenerate the inventory from the project root with:

```sh
python3 scripts/build_data_inventory.py
```

The HTML links assume this directory sits next to `data/`, as in:

```text
scaling_law/
  data/
  shared_data/
```
