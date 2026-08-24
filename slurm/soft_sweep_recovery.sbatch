#!/usr/bin/env bash
#SBATCH --job-name=soft-sweep-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-47
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
CELL_INDEX=$((TASK_INDEX % 12))
REPLICATE_INDEX=$((TASK_INDEX / 12))
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task6_soft_sweep/recovery/task_${TASK_INDEX}"
cd "$PROJECT_ROOT"

SELECTIONS=(0.01 0.02 0.05 0.10 0.01 0.02 0.05 0.10 0.01 0.02 0.05 0.10)
START_AGES=(200 200 200 200 500 500 500 500 1000 1000 1000 1000)
NEUTRAL_DELAYS=(800 1200 2000 3200 1000 1800 3000 4500 1600 2600 4200 6000)
START_FREQUENCIES=(0.01 0.03 0.05 0.10 0.03 0.05 0.10 0.01 0.05 0.10 0.01 0.03)
END_FREQUENCIES=(0.80 0.95 0.99 0.95 0.99 0.80 0.95 0.99 0.95 0.99 0.80 0.95)
MUTATION_GENERATION=$((START_AGES[CELL_INDEX] + NEUTRAL_DELAYS[CELL_INDEX]))
POSITION=$((36500000 + ((CELL_INDEX * 431129 + REPLICATE_INDEX * 71993 + 193733) % 3500000)))

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
  --selection-start-generation "${START_AGES[$CELL_INDEX]}" \
  --minimum-start-frequency "${START_FREQUENCIES[$CELL_INDEX]}" \
  --selection-coefficient "${SELECTIONS[$CELL_INDEX]}" \
  --minimum-end-frequency "${END_FREQUENCIES[$CELL_INDEX]}" \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed $((75226000 + TASK_INDEX))
