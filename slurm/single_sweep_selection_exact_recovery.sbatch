#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-exact-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-31
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
CELL_INDEX=$((TASK_INDEX % 8))
REPLICATE_INDEX=$((TASK_INDEX / 8))
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_exact_recovery/task_${TASK_INDEX}"
cd "$PROJECT_ROOT"

SELECTIONS=(0.02 0.05 0.1 0.2 0.02 0.05 0.1 0.2)
AGES=(4000 1800 900 500 6000 2700 1350 750)
END_FREQUENCIES=(0.95 0.99 0.95 0.99 0.99 0.95 0.99 0.95)
POSITION=$((36500000 + ((CELL_INDEX * 431129 + REPLICATE_INDEX * 71993 + 91733) % 3500000)))

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
  --mutation-generation "${AGES[$CELL_INDEX]}" \
  --selection-coefficient "${SELECTIONS[$CELL_INDEX]}" \
  --minimum-end-frequency "${END_FREQUENCIES[$CELL_INDEX]}" \
  --slim-scaling-factor 1 \
  --slim-burn-in 0.1 \
  --seed $((75225000 + TASK_INDEX))
