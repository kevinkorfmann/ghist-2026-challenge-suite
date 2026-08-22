#!/usr/bin/env bash
#SBATCH --job-name=growth-select-smoke
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/simulate_growth_bgs_selection.py \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --out artifacts/task3_growth_bgs/selection_smoke/fit.json \
  --trees artifacts/task3_growth_bgs/selection_smoke/simulation.trees \
  --slim-path external/SLiM-v4.3/build/slim \
  --chromosome 2L \
  --left 5000000 \
  --right 6000000 \
  --ancestral-size 10000 \
  --generations 5000 \
  --scaling 2.0 \
  --slim-scaling-factor 10 \
  --slim-burn-in 0.1 \
  --samples 42 \
  --seed 7593401 \
  --fit-restarts 32
