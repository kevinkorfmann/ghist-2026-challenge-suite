#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-cal
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="/vast/projects/smathi/cohort/korfmann/singer-amortized"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
SINGER="$REFERENCE_ROOT/build/singer"
CONVERTER="$REFERENCE_ROOT/references/SINGER/SINGER/SINGER/convert_to_tskit"
POLEGON="$REFERENCE_ROOT/references/POLEGON/POLEGON/POLEGON/polegon_master"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/arg_calibration/seed_${SEED:-71001}"
INPUT_PREFIX="$RUN_ROOT/simulated"
SINGER_PREFIX="$RUN_ROOT/singer"
TOPOLOGY_PREFIX="$RUN_ROOT/topology"
POLEGON_INPUT="$RUN_ROOT/polegon_input_8"
POLEGON_OUTPUT="$RUN_ROOT/polegon_8"
mkdir -p "$RUN_ROOT"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/simulate_secondary_contact_arg_calibration.py \
  --prefix "$INPUT_PREFIX" --seed "${SEED:-71001}" \
  --ancestral-size "${ANCESTRAL_SIZE:-18000}" \
  --mainland-size "${MAINLAND_SIZE:-24000}" \
  --island-size "${ISLAND_SIZE:-9000}" \
  --split-time "${SPLIT_TIME:-12000}" \
  --contact-time "${CONTACT_TIME:-1500}" \
  --migration-rate "${MIGRATION_RATE:-8e-5}" \
  --sequence-length "${WINDOW_LENGTH:-100000}"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$INPUT_PREFIX.truth.trees" --out "$RUN_ROOT/truth_hazards.npz"
"$SINGER" \
  -Ne "${ANCESTRAL_SIZE:-18000}" -m 8.4e-8 -r 5.39e-8 \
  -input "$INPUT_PREFIX" -output "$SINGER_PREFIX" \
  -start 0 -end "${WINDOW_LENGTH:-100000}" -ploidy 2 \
  -n 16 -thin 10 -seed "$(( ${SEED:-71001} + 1000 ))"
"$REFERENCE_PYTHON" "$CONVERTER" \
  -input "$SINGER_PREFIX" -output "$TOPOLOGY_PREFIX" -start 8 -end 16
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$RUN_ROOT/topology_*.trees" --out "$RUN_ROOT/singer_hazards.npz"
.venv/bin/python scripts/prepare_polegon_tree_sequence.py \
  --input "${TOPOLOGY_PREFIX}_8.trees" --output "$POLEGON_INPUT.trees"
"$REFERENCE_PYTHON" "$POLEGON" \
  -m 8.4e-8 -input "$POLEGON_INPUT" -output "$POLEGON_OUTPUT" \
  -burn_in 100 -num_samples 100 -thin 10 -scaling_rep 3 \
  -seed "$(( ${SEED:-71001} + 2000 ))" -write_samples
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$POLEGON_OUTPUT.trees" --out "$RUN_ROOT/polegon_hazards.npz"
.venv/bin/python scripts/compare_secondary_contact_arg_features.py \
  --truth "$RUN_ROOT/truth_hazards.npz" \
  --singer "$RUN_ROOT/singer_hazards.npz" \
  --polegon "$RUN_ROOT/polegon_hazards.npz" \
  --out "$RUN_ROOT/comparison.json"
