#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-train-aggregate
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/aggregate_secondary_contact_arg_training.py \
  --root "${ARG_TRAINING_ROOT:-artifacts/task2_secondary_contact/arg_training_v2}" \
  --out data/processed/secondary_contact/arg_training.npz \
  --cells "${CELLS:-2304}"
