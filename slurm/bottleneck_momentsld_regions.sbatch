#!/usr/bin/env bash
#SBATCH --job-name=bn-mld-region
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-49
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python -m ghist_bottleneck.momentsld_reference region \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --block-index "$SLURM_ARRAY_TASK_ID" \
  --out "artifacts/task1_bottleneck/momentsld_regions/region_${SLURM_ARRAY_TASK_ID}.joblib"
