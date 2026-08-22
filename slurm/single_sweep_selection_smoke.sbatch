#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-selection-smoke
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_smoke"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/simulate_single_sweep_selection.py \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --scan "$OUTPUT_ROOT.csv" \
  --slim-path external/SLiM-v4.3/build/slim \
  --length 5000000 \
  --position 2500000 \
  --mutation-generation 1000 \
  --selection-coefficient 0.05 \
  --minimum-end-frequency 0.95 \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed 75221001

