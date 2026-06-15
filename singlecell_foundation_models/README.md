# Single-cell Foundation Models

Reusable workspace for running current single-cell foundation model tooling from the
RMT_sc project.

## Remote Location

Server path:

```sh
/n/holylabs/rongma_lab/Lab/scaling_law/singlecell_foundation_models
```

Project virtual environment:

```sh
/n/holylabs/rongma_lab/Lab/scaling_law/singlecell_foundation_models/.venv
```

On the lab filesystem, `.venv` currently points to the verified environment at:

```sh
/n/home12/jingyuanhu/RMT_sc/singlecell_foundation_models/.venv
```

## Quick Start On The Server

```sh
cd /n/holylabs/rongma_lab/Lab/scaling_law/singlecell_foundation_models
source scripts/activate.sh
python scripts/smoke_test.py
```

To rebuild from scratch:

```sh
RESET_VENV=1 bash scripts/setup_env.sh
```

## What This Sets Up

- `scGPT` from PyPI for embeddings, annotation, integration, perturbation, and
  related workflows.
- `Geneformer` dependencies and a local Hugging Face model checkout under
  `models/geneformer/Geneformer`.
- `UCE` source checkout under `sources/UCE` with the default 4-layer model,
  token table, and protein embeddings under `models/uce/model_files`.
- `scFoundation` source checkout under `sources/scFoundation` with a Hugging
  Face checkpoint mirror under `models/scfoundation/models.ckpt`.
- A general Python single-cell stack: PyTorch, Scanpy, AnnData, Hugging Face
  Transformers/Datasets/Accelerate, Jupyter, and plotting tools.
- scGPT dependency stack including `scvi-tools` and `scib`.

## Directory Layout

```text
singlecell_foundation_models/
  .venv/                 # created on the server, not tracked
  models/                # local checkpoint/model cache
  sources/               # cloned model repositories
  notebooks/             # analysis notebooks
  scripts/
    activate.sh          # load Python module and activate the venv
    setup_env.sh         # create/update the venv
    fetch_models.sh      # clone/download model sources/checkpoints
    smoke_test.py        # import + hardware sanity check
  src/                   # project-specific wrappers/utilities
```

## Notes

- scGPT documents `flash-attn` as optional. It is intentionally not installed by
  default because it is sensitive to the cluster CUDA/compiler stack. Add it only
  inside a GPU job after checking the active CUDA module.
- `torch.cuda.is_available()` is expected to be `False` on the login node. Check
  again inside an allocated GPU job.
- Geneformer’s official Hugging Face instructions recommend `git-lfs` before
  cloning the model repository. The server already had `git-lfs`, and the model
  repository was cloned successfully.
- UCE defaults to the 4-layer checkpoint. The 33-layer checkpoint is larger and
  intentionally not downloaded by default.
- scFoundation inference code calls `.cuda()` in the upstream loader, so run
  full embedding jobs inside an allocated GPU job rather than on the login node.
- Large checkpoints should stay under `models/` on the server, not in git.
