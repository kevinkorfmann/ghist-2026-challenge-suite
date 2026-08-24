#!/usr/bin/env bash
#SBATCH --job-name=bgs-sweeps-localize
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=01:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  VCF="data/raw/sweeps/single_sweep_background_selection/testing/onesweep_bgs.testing.vcf.gz"
  REGION=(105800000 107200000)
  OUT="artifacts/task7_single_sweep_bgs/localization"
else
  VCF="data/raw/sweeps/multiple_sweeps_background_selection/testing/multisweep_bgs.testing.vcf.gz"
  case "$SLURM_ARRAY_TASK_ID" in
    1) REGION=(75500000 77200000) ;;
    2) REGION=(87400000 89100000) ;;
    3) REGION=(89600000 91400000) ;;
    4) REGION=(129100000 130800000) ;;
    5) REGION=(137700000 139300000) ;;
    6) REGION=(154900000 156600000) ;;
  esac
  OUT="artifacts/task8_multiple_sweeps_bgs/localization/region_${SLURM_ARRAY_TASK_ID}"
fi

PYTHONPATH=src .venv/bin/python scripts/localize_single_sweep.py \
  --vcf "$VCF" \
  --region "${REGION[@]}" \
  --grid-step 1000 \
  --out-dir "$OUT" \
  --summary "$OUT/summary.json"
