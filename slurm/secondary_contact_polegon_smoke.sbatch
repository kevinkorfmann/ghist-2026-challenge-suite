#!/usr/bin/env bash
#SBATCH --job-name=sc-polegon
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="/vast/projects/smathi/cohort/korfmann/singer-amortized"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
POLEGON="$REFERENCE_ROOT/references/POLEGON/POLEGON/POLEGON/polegon_master"
WINDOW_INDEX="${WINDOW_INDEX:-500}"
DRAW="${DRAW:-8}"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/singer/window_${WINDOW_INDEX}"
INPUT_PREFIX="$RUN_ROOT/topology_${DRAW}"
OUTPUT_PREFIX="$RUN_ROOT/polegon_${DRAW}"
cd "$PROJECT_ROOT"
test -s "$INPUT_PREFIX.trees"
"$REFERENCE_PYTHON" "$POLEGON" \
  -m 8.4e-8 \
  -input "$INPUT_PREFIX" -output "$OUTPUT_PREFIX" \
  -burn_in "${BURNIN:-100}" \
  -num_samples "${SAMPLES:-100}" \
  -thin "${THIN:-10}" \
  -scaling_rep "${SCALING_REPS:-3}" \
  -seed "${SEED:-68001}" -write_samples
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$OUTPUT_PREFIX.trees" \
  --out "$RUN_ROOT/polegon_${DRAW}_coalescence_hazards.npz"
