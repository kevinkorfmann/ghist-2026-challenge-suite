#!/usr/bin/env bash
#SBATCH --job-name=sc-singer-tile
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
SINGER="$REFERENCE_ROOT/build/singer"
CONVERTER="$REFERENCE_ROOT/references/SINGER/SINGER/SINGER/convert_to_tskit"
REFERENCE_PYTHON="$REFERENCE_ROOT/.venv-cluster/bin/python"
printf -v TILE "%04d" "$SLURM_ARRAY_TASK_ID"
INPUT_ROOT="$PROJECT_ROOT/data/processed/secondary_contact/singer_tiles_100kb/tile_${TILE}"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/singer_tiles_100kb/tile_${TILE}"
INPUT_PREFIX="$INPUT_ROOT/observed"
OUTPUT_PREFIX="$RUN_ROOT/singer"
mkdir -p "$RUN_ROOT"
cd "$PROJECT_ROOT"
test -s "$INPUT_PREFIX.vcf"
"$SINGER" \
  -Ne "${EFFECTIVE_SIZE:-18000}" -m 8.4e-8 -r 5.39e-8 \
  -input "$INPUT_PREFIX" -output "$OUTPUT_PREFIX" \
  -start 0 -end 100000 -ploidy 2 \
  -n "${DRAWS:-12}" -thin "${THIN:-10}" \
  -seed "$(( ${SEED_BASE:-75000} + SLURM_ARRAY_TASK_ID ))"
"$REFERENCE_PYTHON" "$CONVERTER" \
  -input "$OUTPUT_PREFIX" -output "$RUN_ROOT/topology" \
  -start "${BURNIN_DRAWS:-8}" -end "${DRAWS:-12}"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$RUN_ROOT/topology_*.trees" \
  --out "$RUN_ROOT/coalescence_hazards.npz"
