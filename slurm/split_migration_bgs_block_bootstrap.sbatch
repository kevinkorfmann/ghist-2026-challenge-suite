#!/usr/bin/env bash
#SBATCH --job-name=splitmig-block-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-63%16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/bootstrap_split_migration_bgs_blocks.py \
  --blocks artifacts/task4_split_migration_bgs/blocks_1mb.npz \
  --initial-fit experiments/split_migration_bgs/baseline_p50_40.json \
  --replicate "$SLURM_ARRAY_TASK_ID" \
  --restarts 8 \
  --out "artifacts/task4_split_migration_bgs/bootstrap_1mb/task_${SLURM_ARRAY_TASK_ID}.json"

