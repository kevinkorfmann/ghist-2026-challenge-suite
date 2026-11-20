#!/usr/bin/env bash
#SBATCH --job-name=bgs-sweeps-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  VCF="data/raw/sweeps/single_sweep_background_selection/testing/onesweep_bgs.testing.vcf.gz"
  OUT="artifacts/task7_single_sweep_bgs"
else
  VCF="data/raw/sweeps/multiple_sweeps_background_selection/testing/multisweep_bgs.testing.vcf.gz"
  OUT="artifacts/task8_multiple_sweeps_bgs"
fi

PYTHONPATH=src .venv/bin/python scripts/audit_single_sweep.py \
  --vcf "$VCF" \
  --accessible-bed data/raw/shared/sweeps/testing.accessible.bed \
  --score-field composite_score \
  --out-dir "$OUT/baseline_scans" \
  --summary "$OUT/audit_summary.json"
