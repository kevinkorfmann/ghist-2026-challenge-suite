#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-selection-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-63
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_recovery/task_${TASK_INDEX}"
cd "$PROJECT_ROOT"

SELECTIONS=(0.02 0.05 0.1 0.2)
BASE_AGES=(4000 1800 900 500)
END_FREQUENCIES=(0.9 0.95 0.99 0.999)
AGE_NUMERATORS=(3 4 6 8)
SELECTION_INDEX=$((TASK_INDEX % 4))
AGE_INDEX=$(((TASK_INDEX / 4) % 4))
END_INDEX=$(((TASK_INDEX / 16) % 4))
POSITION=$((1000000 + ((TASK_INDEX * 7919 + 104729) % 3000000)))
MUTATION_GENERATION=$((BASE_AGES[SELECTION_INDEX] * AGE_NUMERATORS[AGE_INDEX] / 4))

PYTHONPATH=src .venv/bin/python scripts/simulate_single_sweep_selection.py \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --scan "$OUTPUT_ROOT.csv" \
  --slim-path external/SLiM-v4.3/build/slim \
  --length 5000000 \
  --position "$POSITION" \
  --mutation-generation "$MUTATION_GENERATION" \
  --selection-coefficient "${SELECTIONS[$SELECTION_INDEX]}" \
  --minimum-end-frequency "${END_FREQUENCIES[$END_INDEX]}" \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed $((75222000 + TASK_INDEX))

