#!/usr/bin/env bash
#SBATCH --job-name=splitmig-bgs-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PROJECTIONS=("94 68" "70 50" "50 40")
read -r P1 P2 <<< "${PROJECTIONS[$SLURM_ARRAY_TASK_ID]}"

PYTHONPATH=src .venv/bin/python scripts/fit_split_migration_bgs_baseline.py \
  --vcf data/raw/demography/split_migration_background_selection/testing/split_migration.testing.vcf.gz \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --projection "$P1" "$P2" \
  --restarts 64 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --out "artifacts/task4_split_migration_bgs/baseline_p${P1}_${P2}.json"

