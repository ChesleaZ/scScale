#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

mkdir -p sources models/scgpt models/geneformer models/scfoundation models/huggingface

if [ ! -d sources/scGPT ]; then
  git clone https://github.com/bowang-lab/scGPT.git sources/scGPT
else
  git -C sources/scGPT pull --ff-only
fi

if [ ! -d sources/scFoundation ]; then
  git clone https://github.com/biomap-research/scFoundation.git sources/scFoundation
else
  git -C sources/scFoundation pull --ff-only
fi

if command -v git-lfs >/dev/null 2>&1; then
  git lfs install
  if [ ! -d models/geneformer/Geneformer ]; then
    git clone https://huggingface.co/ctheodoris/Geneformer models/geneformer/Geneformer
  else
    git -C models/geneformer/Geneformer pull --ff-only
  fi
else
  echo "git-lfs not found; skipping full Geneformer repository clone."
  echo "Use huggingface_hub or load ctheodoris/Geneformer from Python when needed."
fi

cat <<'MSG'

Model/source fetch complete.

scGPT pretrained checkpoints are distributed through the links in the official
scGPT model zoo. Download chosen checkpoint folders into:
  models/scgpt/

scFoundation weights and examples should go into:
  models/scfoundation/

MSG

