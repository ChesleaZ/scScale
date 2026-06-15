#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if command -v module >/dev/null 2>&1; then
  module load python/3.10.13-fasrc01
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load python/3.10.13-fasrc01
fi

source "${ROOT}/.venv/bin/activate"
export HF_HOME="${ROOT}/models/huggingface"
export HF_DATASETS_CACHE="${HF_HOME}/datasets"
export SCFM_HOME="${ROOT}"

echo "Activated single-cell foundation model environment:"
echo "  ${VIRTUAL_ENV}"
python --version
