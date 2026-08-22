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
SEED=$((2026082100 + SLURM_ARRAY_TASK_ID))
.venv/bin/python -m ghist_bottleneck.cli simulate \
  --observed data/processed/bottleneck/observed.blocks.npz \
  --out "data/simulations/bottleneck/part_${SLURM_ARRAY_TASK_ID}.npz" \
  --n-sims 1000 \
  --seed "$SEED"
