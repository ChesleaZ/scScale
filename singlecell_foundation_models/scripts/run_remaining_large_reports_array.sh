#!/usr/bin/env bash
#SBATCH --job-name=scfm_remaining_reports
#SBATCH --account=rongma_lab
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-5%2
#SBATCH --output=/n/holylabs/rongma_lab/Lab/scaling_law/reports/logs/remaining_report_%A_%a.out
#SBATCH --error=/n/holylabs/rongma_lab/Lab/scaling_law/reports/logs/remaining_report_%A_%a.err

set -euo pipefail

ROOT="${ROOT:-/n/holylabs/rongma_lab/Lab/scaling_law}"
SCRIPT="${ROOT}/singlecell_foundation_models/scripts/dataset_qc_pca_report.py"
DATA="${ROOT}/data"
REPORTS="${ROOT}/reports/datasets"

mkdir -p "${REPORTS}" "${ROOT}/reports/logs"
source "${ROOT}/singlecell_foundation_models/scripts/activate.sh"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

run_report() {
  local dataset_name="$1"
  local slug="$2"
  local input_type="$3"
  local input_path="$4"
  local technology="$5"
  local assay="$6"
  shift 6

  echo "[$(date)] Starting ${dataset_name}"
  python "${SCRIPT}" \
    --dataset-name "${dataset_name}" \
    --input-type "${input_type}" \
    --input-path "${input_path}" \
    --out-dir "${REPORTS}/${slug}" \
    --technology "${technology}" \
    --assay "${assay}" \
    --max-cells 20000 \
  "$@"
  echo "[$(date)] Finished ${dataset_name}"
}

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0)
    run_report "Norman-Weissman 2019" "NormanWeissman2019" "h5ad" \
      "${DATA}/NormanWeissman2019/NormanWeissman2019_filtered.h5ad" \
      "Perturb-seq" "Gene expression"
    ;;
  1)
    run_report "Replogle RPE1" "Replogle_RPE1" "h5ad" \
      "${DATA}/Replogle_RPE1/rpe1_raw_singlecell_01.h5ad" \
      "Perturb-seq" "Gene expression"
    ;;
  2)
    run_report "Jurkat" "Jurkat" "10x_mtx" \
      "${DATA}/Jurkat/sample_filtered_feature_bc_matrix" \
      "10x Genomics" "Gene expression"
    ;;
  3)
    run_report "K562" "K562" "10x_mtx" \
      "${DATA}/K562/sample_filtered_feature_bc_matrix" \
      "10x Genomics" "Gene expression"
    ;;
  4)
    run_report "Human Lung Cell Atlas core" "HLCA_core" "h5ad" \
      "${DATA}/HLCA/HLCA_core_CELLxGENE_2025-11-08.h5ad" \
      "cellxgene h5ad" "Gene expression"
    ;;
  5)
    run_report "Parse 1M PBMC T1D" "Parse_1M_PBMC_T1D" "parse_zip" \
      "${DATA}/Parse_1M_PBMC_T1D/1M_PBMC_T1D_Parse.zip" \
      "Parse Biosciences" "Gene expression"
    ;;
  *)
    echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac
