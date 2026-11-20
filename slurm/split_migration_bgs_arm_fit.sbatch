#!/usr/bin/env bash
#SBATCH --job-name=splitmig-arm-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

ARMS=(2L 2R 3L 3R)
ARM_INDEX=$((SLURM_ARRAY_TASK_ID % 4))
ARM="${ARMS[$ARM_INDEX]}"
if (( SLURM_ARRAY_TASK_ID < 4 )); then
  MODE=arm_only
else
  MODE=leave_one_out
fi

PYTHONPATH=src .venv/bin/python scripts/fit_split_migration_bgs_arm.py \
  --blocks artifacts/task4_split_migration_bgs/blocks_1mb.npz \
  --initial-fit experiments/split_migration_bgs/baseline_p50_40.json \
  --arm "$ARM" \
  --mode "$MODE" \
  --restarts 32 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --out "artifacts/task4_split_migration_bgs/arm_${MODE}_${ARM}.json"
