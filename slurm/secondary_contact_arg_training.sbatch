#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-train-data
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="/vast/projects/smathi/cohort/korfmann/singer-amortized"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
SINGER="$REFERENCE_ROOT/build/singer"
CONVERTER="$REFERENCE_ROOT/references/SINGER/SINGER/SINGER/convert_to_tskit"
POLEGON="$REFERENCE_ROOT/references/POLEGON/POLEGON/POLEGON/polegon_master"
printf -v CELL "%05d" "$SLURM_ARRAY_TASK_ID"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/arg_training/cell_${CELL}"
INPUT_PREFIX="$RUN_ROOT/simulated"
SINGER_PREFIX="$RUN_ROOT/singer"
TOPOLOGY_PREFIX="$RUN_ROOT/topology"
POLEGON_INPUT="$RUN_ROOT/polegon_input_8"
POLEGON_OUTPUT="$RUN_ROOT/polegon_8"
mkdir -p "$RUN_ROOT"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/simulate_secondary_contact_arg_design.py \
  --prefix "$INPUT_PREFIX" \
  --cell-id "$SLURM_ARRAY_TASK_ID" \
  --points "${POINTS:-768}" \
  --replicates "${REPLICATES:-3}" \
  --center \
  "${ANCESTRAL:-16585.6}" "${MAINLAND:-18883.2}" "${ISLAND:-51091.5}" \
  "${SPLIT:-4622.4}" "${CONTACT:-1188.0}" "${MIGRATION:-0.00024268}" \
  --seed "${DESIGN_SEED:-20260903}"
"$SINGER" \
  -Ne "${SINGER_SIZE:-18000}" -m 8.4e-8 -r 5.39e-8 \
  -input "$INPUT_PREFIX" -output "$SINGER_PREFIX" \
  -start 0 -end 100000 -ploidy 2 \
  -n 12 -thin 10 -seed "$(( ${SINGER_SEED_BASE:-91000} + SLURM_ARRAY_TASK_ID ))"
"$REFERENCE_PYTHON" "$CONVERTER" \
  -input "$SINGER_PREFIX" -output "$TOPOLOGY_PREFIX" -start 8 -end 12
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$RUN_ROOT/topology_*.trees" --out "$RUN_ROOT/singer_hazards.npz"
.venv/bin/python scripts/prepare_polegon_tree_sequence.py \
  --input "${TOPOLOGY_PREFIX}_8.trees" --output "$POLEGON_INPUT.trees"
"$REFERENCE_PYTHON" "$POLEGON" \
  -m 8.4e-8 -input "$POLEGON_INPUT" -output "$POLEGON_OUTPUT" \
  -burn_in 100 -num_samples 100 -thin 10 -scaling_rep 3 \
  -seed "$(( ${POLEGON_SEED_BASE:-93000} + SLURM_ARRAY_TASK_ID ))"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$POLEGON_OUTPUT.trees" --out "$RUN_ROOT/polegon_hazards.npz"
