#!/usr/bin/env bash
#SBATCH --job-name=sc-sfs
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
POINTS="${POINTS:-2000}"
POINTS_PER_TASK="${POINTS_PER_TASK:-20}"
POINT_START=$((SLURM_ARRAY_TASK_ID * POINTS_PER_TASK))
DESIGN="${DESIGN:-broad}"
OUTPUT_ROOT="${OUTPUT_ROOT:-data/processed/secondary_contact/sfs_chunks}"
EXTRA_ARGS=()
if [[ "$DESIGN" == "focused" ]]; then
  EXTRA_ARGS+=(
    --center
    "${ANCESTRAL:-16585}" "${MAINLAND:-18880}" "${ISLAND:-51090}"
    "${SPLIT:-4622}" "${CONTACT:-1188}" "${MIGRATION:-0.000243}"
    --point-id-offset "${POINT_ID_OFFSET:-1000000}"
  )
fi
.venv/bin/python scripts/generate_secondary_contact_sfs_dataset.py \
  --out "$OUTPUT_ROOT/chunk_${SLURM_ARRAY_TASK_ID}.npz" \
  --point-start "$POINT_START" \
  --point-count "$POINTS_PER_TASK" \
  --points "$POINTS" \
  --replicates "${REPLICATES:-6}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --design "$DESIGN" \
  "${EXTRA_ARGS[@]}"
