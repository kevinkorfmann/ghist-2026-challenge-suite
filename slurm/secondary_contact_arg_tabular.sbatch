#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-tabular
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/train_secondary_contact_arg_tabular.py \
  --dataset data/processed/secondary_contact/arg_training.npz \
  --out artifacts/task2_secondary_contact/model_arg_tabular_core.joblib \
  --domain-filter 0 \
  --workers "$SLURM_CPUS_PER_TASK"
.venv/bin/python scripts/predict_secondary_contact_arg_tabular.py \
  --features artifacts/task2_secondary_contact/tiled_arg_features.npz \
  --model artifacts/task2_secondary_contact/model_arg_tabular_core.joblib \
  --out artifacts/task2_secondary_contact/arg_tabular_core_testing_prediction.json \
  --tiles-per-block 20
