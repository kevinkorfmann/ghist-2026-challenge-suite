#!/usr/bin/env bash
#SBATCH --job-name=bn-chrom-val
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-19
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
CELLS_PER_TASK=10
CELL_START=$((SLURM_ARRAY_TASK_ID * CELLS_PER_TASK))
SEED=$((2026100100 + SLURM_ARRAY_TASK_ID))
.venv/bin/python -m ghist_bottleneck.cli simulate-chromosomes \
  --observed data/processed/task1_bottleneck/observed.blocks.npz \
  --out "data/simulations/bottleneck/chromosome_${SLURM_ARRAY_TASK_ID}.npz" \
  --cells "$CELLS_PER_TASK" \
  --blocks-per-cell 50 \
  --cell-start "$CELL_START" \
  --seed "$SEED"
