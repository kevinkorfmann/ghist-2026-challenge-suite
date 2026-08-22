#!/usr/bin/env bash
#SBATCH --job-name=splitmig-select-grid
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-15%8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

ANCESTRAL=(10000 12000 15000 18000)
AFR=(12000 24000 15000 30000)
EUR=(18000 16000 30000 22000)
GENERATIONS=(2500 8000 3500 6500)
MIGRATION=(0.00008 0.00025 0.00012 0.00020)
CHROMOSOMES=(2L 2R 3L 3R)
LEFTS=(5000000 10000000 10000000 15000000)

ARM_INDEX=$((SLURM_ARRAY_TASK_ID % 4))
POINT_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
CHROMOSOME="${CHROMOSOMES[$ARM_INDEX]}"
LEFT="${LEFTS[$ARM_INDEX]}"
RIGHT=$((LEFT + 1000000))
OUTPUT_ROOT="artifacts/task4_split_migration_bgs/selection_grid/point_${POINT_INDEX}/task_${ARM_INDEX}"

PYTHONPATH=src .venv/bin/python scripts/simulate_split_migration_bgs_selection.py \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome "$CHROMOSOME" \
  --left "$LEFT" \
  --right "$RIGHT" \
  --ancestral-size "${ANCESTRAL[$POINT_INDEX]}" \
  --afr-size "${AFR[$POINT_INDEX]}" \
  --eur-size "${EUR[$POINT_INDEX]}" \
  --generations "${GENERATIONS[$POINT_INDEX]}" \
  --migration-rate "${MIGRATION[$POINT_INDEX]}" \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed $((75196000 + SLURM_ARRAY_TASK_ID)) \
  --skip-fit

