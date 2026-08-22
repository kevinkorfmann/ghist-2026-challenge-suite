#!/usr/bin/env bash
#SBATCH --job-name=splitmig-distance-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

STRATA=(coding near far all_noncoding)
STRATUM="${STRATA[$SLURM_ARRAY_TASK_ID]}"

PYTHONPATH=src .venv/bin/python scripts/fit_split_migration_bgs_baseline.py \
  --vcf data/raw/demography/split_migration_background_selection/testing/split_migration.testing.vcf.gz \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --projection 50 40 \
  --stratum "$STRATUM" \
  --restarts 32 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --out "artifacts/task4_split_migration_bgs/distance_${STRATUM}_p50_40.json"
