#!/usr/bin/env bash
#SBATCH --job-name=growth-select-compare
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/compare_growth_bgs_selection_smoke.py \
  --simulation artifacts/task3_growth_bgs/selection_smoke/fit.json \
  --sfs artifacts/task3_growth_bgs/selection_smoke/fit.sfs.npz \
  --out artifacts/task3_growth_bgs/selection_smoke/convention_comparison.json \
  --restarts 32
