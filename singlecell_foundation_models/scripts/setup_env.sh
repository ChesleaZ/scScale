#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

if command -v module >/dev/null 2>&1; then
  module load python/3.10.13-fasrc01
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load python/3.10.13-fasrc01
fi

if [ "${RESET_VENV:-0}" = "1" ]; then
  rm -rf .venv
fi

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip "setuptools==80.9.0" wheel
python -m pip install -r requirements.txt
python -m ipykernel install --user --name rmt-sc-fm --display-name "RMT sc foundation models"

mkdir -p models/huggingface models/scgpt models/geneformer models/scfoundation sources notebooks src

python scripts/smoke_test.py
