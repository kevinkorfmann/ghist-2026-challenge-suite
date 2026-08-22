#!/usr/bin/env bash
#SBATCH --job-name=sc-polegon-tile
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="/vast/projects/smathi/cohort/korfmann/singer-amortized"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
POLEGON="$REFERENCE_ROOT/references/POLEGON/POLEGON/POLEGON/polegon_master"
printf -v TILE "%04d" "$SLURM_ARRAY_TASK_ID"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/singer_tiles_100kb/tile_${TILE}"
TOPOLOGY="$RUN_ROOT/topology_8.trees"
INPUT_PREFIX="$RUN_ROOT/polegon_input_8"
OUTPUT_PREFIX="$RUN_ROOT/polegon_8"
cd "$PROJECT_ROOT"
test -s "$TOPOLOGY"
.venv/bin/python scripts/prepare_polegon_tree_sequence.py \
  --input "$TOPOLOGY" --output "$INPUT_PREFIX.trees"
"$REFERENCE_PYTHON" "$POLEGON" \
  -m 8.4e-8 -input "$INPUT_PREFIX" -output "$OUTPUT_PREFIX" \
  -burn_in 100 -num_samples 100 -thin 10 -scaling_rep 3 \
  -seed "$(( ${SEED_BASE:-85000} + SLURM_ARRAY_TASK_ID ))"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$OUTPUT_PREFIX.trees" \
  --out "$RUN_ROOT/polegon_coalescence_hazards.npz"
