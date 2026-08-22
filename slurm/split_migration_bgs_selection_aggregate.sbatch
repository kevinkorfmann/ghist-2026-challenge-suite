#!/usr/bin/env bash
#SBATCH --job-name=splitmig-select-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/aggregate_split_migration_bgs_selection.py \
  --input-dir artifacts/task4_split_migration_bgs/selection_candidate \
  --out artifacts/task4_split_migration_bgs/selection_candidate_aggregate.json \
  --restarts 64 \
  --workers "$SLURM_CPUS_PER_TASK"
