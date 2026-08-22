#!/usr/bin/env bash
#SBATCH --job-name=splitmig-select-candidate
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

CHROMOSOMES=(2L 2R 3L 3R)
LEFTS=(5000000 10000000 10000000 15000000)
CHROMOSOME="${CHROMOSOMES[$SLURM_ARRAY_TASK_ID]}"
LEFT="${LEFTS[$SLURM_ARRAY_TASK_ID]}"
RIGHT=$((LEFT + 1000000))
OUTPUT_ROOT="artifacts/task4_split_migration_bgs/selection_candidate/task_${SLURM_ARRAY_TASK_ID}"

PYTHONPATH=src .venv/bin/python scripts/simulate_split_migration_bgs_selection.py \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --out "$OUTPUT_ROOT.json" \
  --trees "$OUTPUT_ROOT.trees" \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome "$CHROMOSOME" \
  --left "$LEFT" \
  --right "$RIGHT" \
  --ancestral-size 15000 \
  --afr-size 20000 \
  --eur-size 30000 \
  --generations 5000 \
  --migration-rate 0.00015 \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --seed $((75194000 + SLURM_ARRAY_TASK_ID)) \
  --fit-restarts 16

