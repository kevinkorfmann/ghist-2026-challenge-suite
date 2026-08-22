#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-chr7-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_chr7_recovery/task_${TASK_INDEX}"
cd "$PROJECT_ROOT"

SELECTIONS=(0.02 0.05 0.1 0.2)
BASE_AGES=(4000 1800 900 500)
AGE_NUMERATORS=(3 4 6 8)
SELECTION_INDEX=$((TASK_INDEX % 4))
AGE_INDEX=$((TASK_INDEX / 4))
POSITION=$((36500000 + ((TASK_INDEX * 104729 + 65537) % 3500000)))
MUTATION_GENERATION=$((BASE_AGES[SELECTION_INDEX] * AGE_NUMERATORS[AGE_INDEX] / 4))
if (( TASK_INDEX % 2 == 0 )); then
  END_FREQUENCY=0.95
else
  END_FREQUENCY=0.99
fi

PYTHONPATH=src .venv/bin/python scripts/simulate_single_sweep_selection.py \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --scan "$OUTPUT_ROOT.csv" \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome 7 \
  --left 35750000 \
  --genetic-map HapMapII_GRCh38 \
  --length 5000000 \
  --position "$POSITION" \
  --mutation-generation "$MUTATION_GENERATION" \
  --selection-coefficient "${SELECTIONS[$SELECTION_INDEX]}" \
  --minimum-end-frequency "$END_FREQUENCY" \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed $((75223000 + TASK_INDEX))

