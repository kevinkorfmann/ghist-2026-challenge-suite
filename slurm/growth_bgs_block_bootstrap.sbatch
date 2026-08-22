#!/usr/bin/env bash
#SBATCH --job-name=growth-bgs-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-31%16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

REPLICATES_PER_TASK="${REPLICATES_PER_TASK:-4}"
START=$((SLURM_ARRAY_TASK_ID * REPLICATES_PER_TASK))
PYTHONPATH=src .venv/bin/python scripts/bootstrap_growth_bgs_blocks.py \
  --blocks artifacts/task3_growth_bgs/middle_distance_blocks_5mb.npz \
  --out "artifacts/task3_growth_bgs/block_bootstrap/task_${SLURM_ARRAY_TASK_ID}.json" \
  --replicate-start "$START" \
  --replicates "$REPLICATES_PER_TASK" \
  --projection 84 \
  --restarts 12 \
  --maxiter 400
