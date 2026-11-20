#!/usr/bin/env bash
#SBATCH --job-name=growth-select-stress-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/fit_growth_bgs_selection_stress.py \
  --stress-dir artifacts/task3_growth_bgs/selection_stress \
  --baseline-dir artifacts/task3_growth_bgs/selection_grid \
  --out artifacts/task3_growth_bgs/selection_stress_aggregate.json \
  --restarts 64 \
  --workers "$SLURM_CPUS_PER_TASK"
