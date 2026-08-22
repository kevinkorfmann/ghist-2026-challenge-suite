#!/usr/bin/env bash
#SBATCH --job-name=bn-multi-sfs
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-89
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
CELLS_PER_TASK="${CELLS_PER_TASK:-10}"
TASK_OFFSET="${TASK_OFFSET:-0}"
OUT_PREFIX="${OUT_PREFIX:-multiscale_sfs}"
CELL_START=$(((SLURM_ARRAY_TASK_ID + TASK_OFFSET) * CELLS_PER_TASK))
.venv/bin/python scripts/run_multiscale_sfs_recovery.py \
  --out "data/simulations/bottleneck/multiscale_precision/${OUT_PREFIX}_${SLURM_ARRAY_TASK_ID}.npz" \
  --cell-start "$CELL_START" \
  --cells "$CELLS_PER_TASK" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --blocks 50 \
  --restarts 8 \
  --seed 20261117
