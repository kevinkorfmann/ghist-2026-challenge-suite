#!/usr/bin/env bash
#SBATCH --job-name=sc-sim
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
POINTS="${POINTS:-600}"
REPLICATES="${REPLICATES:-6}"
CELLS_PER_TASK="${CELLS_PER_TASK:-8}"
CELL_START=$((SLURM_ARRAY_TASK_ID * CELLS_PER_TASK))
.venv/bin/python scripts/generate_secondary_contact_dataset.py \
  --out "data/processed/secondary_contact/chunks/chunk_${SLURM_ARRAY_TASK_ID}.npz" \
  --cell-start "$CELL_START" \
  --cells "$CELLS_PER_TASK" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --points "$POINTS" \
  --replicates "$REPLICATES" \
  --sequence-length "${SEQUENCE_LENGTH:-100000000}"
