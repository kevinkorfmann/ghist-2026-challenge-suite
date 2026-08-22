#!/usr/bin/env bash
#SBATCH --job-name=bn-sim
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-99
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
TASK_OFFSET="${TASK_OFFSET:-0}"
SEED_BASE="${SEED_BASE:-2026082100}"
SHARD_PREFIX="${SHARD_PREFIX:-part}"
SIMS_PER_TASK="${SIMS_PER_TASK:-1000}"
TASK_INDEX=$((TASK_OFFSET + SLURM_ARRAY_TASK_ID))
SEED=$((SEED_BASE + TASK_INDEX))
.venv/bin/python -m ghist_bottleneck.cli simulate \
  --observed data/processed/task1_bottleneck/observed.blocks.npz \
  --out "data/simulations/bottleneck/${SHARD_PREFIX}_${TASK_INDEX}.npz" \
  --n-sims "$SIMS_PER_TASK" \
  --seed "$SEED"
