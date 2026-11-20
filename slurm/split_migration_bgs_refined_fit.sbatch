#!/usr/bin/env bash
#SBATCH --job-name=splitmig-refined-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/fit_split_migration_bgs_baseline.py \
  --vcf data/raw/demography/split_migration_background_selection/testing/split_migration.testing.vcf.gz \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --projection 50 40 \
  --stratum middle \
  --restarts 32 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --integration-dt-fac 0.005 \
  --adaptive-integration \
  --initial-fit experiments/split_migration_bgs/baseline_p50_40.json \
  --out artifacts/task4_split_migration_bgs/refined_middle_p50_40.json
