#!/usr/bin/env bash
#SBATCH --job-name=growth-select-grid
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-35%8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

GENERATIONS=(2500 5000 10000)
SCALINGS=(1.5 2.0 3.0)
CHROMOSOMES=(2L 2R 3L 3R)
LEFTS=(5000000 10000000 10000000 15000000)
WINDOW_INDEX=$((SLURM_ARRAY_TASK_ID % 4))
PARAMETER_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
GENERATION_INDEX=$((PARAMETER_INDEX / 3))
SCALING_INDEX=$((PARAMETER_INDEX % 3))
CHROMOSOME="${CHROMOSOMES[$WINDOW_INDEX]}"
LEFT="${LEFTS[$WINDOW_INDEX]}"
RIGHT=$((LEFT + 1000000))
GENERATIONS_VALUE="${GENERATIONS[$GENERATION_INDEX]}"
SCALING_VALUE="${SCALINGS[$SCALING_INDEX]}"
OUTPUT_ROOT="artifacts/task3_growth_bgs/selection_grid/task_${SLURM_ARRAY_TASK_ID}"

PYTHONPATH=src .venv/bin/python scripts/simulate_growth_bgs_selection.py \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome "$CHROMOSOME" \
  --left "$LEFT" \
  --right "$RIGHT" \
  --ancestral-size 10000 \
  --generations "$GENERATIONS_VALUE" \
  --scaling "$SCALING_VALUE" \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --samples 42 \
  --seed $((7593500 + SLURM_ARRAY_TASK_ID)) \
  --fit-restarts 20
