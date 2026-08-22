#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-predict
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/predict_secondary_contact_arg.py \
  --features artifacts/task2_secondary_contact/tiled_arg_features.npz \
  --checkpoints 'artifacts/task2_secondary_contact/model_arg_*_*.pt' \
  --out artifacts/task2_secondary_contact/arg_ensemble_testing_prediction.json \
  --tiles-per-block 20 \
  --device cuda
