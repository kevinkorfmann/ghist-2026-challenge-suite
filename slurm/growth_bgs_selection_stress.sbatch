#!/usr/bin/env bash
#SBATCH --job-name=growth-select-stress
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-11%6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

CHROMOSOMES=(2L 2R 3L 3R)
LEFTS=(5000000 10000000 10000000 15000000)
SCALING_FACTORS=(5 20 10)
BURN_INS=(0.1 0.1 0.5)
WINDOW_INDEX=$((SLURM_ARRAY_TASK_ID % 4))
CONDITION_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
CHROMOSOME="${CHROMOSOMES[$WINDOW_INDEX]}"
LEFT="${LEFTS[$WINDOW_INDEX]}"
RIGHT=$((LEFT + 1000000))
SIMULATION_SCALE="${SCALING_FACTORS[$CONDITION_INDEX]}"
BURN_IN="${BURN_INS[$CONDITION_INDEX]}"
OUTPUT_ROOT="artifacts/task3_growth_bgs/selection_stress/task_${SLURM_ARRAY_TASK_ID}"

PYTHONPATH=src .venv/bin/python scripts/simulate_growth_bgs_selection.py \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome "$CHROMOSOME" \
  --left "$LEFT" \
  --right "$RIGHT" \
  --ancestral-size 10000 \
  --generations 5000 \
  --scaling 2.0 \
  --slim-scaling-factor "$SIMULATION_SCALE" \
  --slim-burn-in "$BURN_IN" \
  --samples 42 \
  --seed $((7595000 + SLURM_ARRAY_TASK_ID)) \
  --fit-restarts 20
