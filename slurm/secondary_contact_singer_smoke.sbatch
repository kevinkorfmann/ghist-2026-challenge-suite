#!/usr/bin/env bash
#SBATCH --job-name=sc-singer
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
SINGER="$REFERENCE_ROOT/build/singer"
CONVERTER="$REFERENCE_ROOT/references/SINGER/SINGER/SINGER/convert_to_tskit"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
WINDOW_INDEX="${WINDOW_INDEX:-50}"
WINDOW_LENGTH="${WINDOW_LENGTH:-100000}"
START=$((WINDOW_INDEX * WINDOW_LENGTH))
END=$((START + WINDOW_LENGTH))
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/singer/window_${WINDOW_INDEX}"
INPUT_PREFIX="$RUN_ROOT/observed"
OUTPUT_PREFIX="$RUN_ROOT/singer"
mkdir -p "$RUN_ROOT"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/prepare_secondary_contact_singer_window.py \
  --vcf data/raw/demography/secondary_contact/testing/secondary_contact.testing.vcf.gz \
  --out "$INPUT_PREFIX.vcf" --start "$START" --end "$END"
test -x "$SINGER"
test -f "$CONVERTER"
"$SINGER" \
  -Ne "${EFFECTIVE_SIZE:-10000}" \
  -m 8.4e-8 -r 5.39e-8 \
  -input "$INPUT_PREFIX" -output "$OUTPUT_PREFIX" \
  -start 0 -end "$WINDOW_LENGTH" -ploidy 2 \
  -n "${DRAWS:-16}" -thin "${THIN:-10}" -seed "${SEED:-65001}"
"$REFERENCE_PYTHON" "$CONVERTER" \
  -input "$OUTPUT_PREFIX" -output "$RUN_ROOT/topology" \
  -start "${BURNIN_DRAWS:-8}" -end "${DRAWS:-16}"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$RUN_ROOT/topology_*.trees" \
  --out "$RUN_ROOT/coalescence_hazards.npz"
