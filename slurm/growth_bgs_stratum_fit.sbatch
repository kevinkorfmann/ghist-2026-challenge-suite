#!/usr/bin/env bash
#SBATCH --job-name=growth-bgs-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PROJECTION="${PROJECTION:-84}"
RESTARTS="${RESTARTS:-64}"
MAXITER="${MAXITER:-600}"
OUTPUT="${OUTPUT:-artifacts/task3_growth_bgs/stratum_fit_p${PROJECTION}.json}"

PYTHONPATH=src .venv/bin/python scripts/fit_growth_bgs_strata.py \
  --audit artifacts/task3_growth_bgs/observed_resource_audit.json \
  --out "$OUTPUT" \
  --projection "$PROJECTION" \
  --restarts "$RESTARTS" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --maxiter "$MAXITER"
