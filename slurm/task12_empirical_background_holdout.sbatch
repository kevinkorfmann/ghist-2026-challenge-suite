#!/usr/bin/env bash
#SBATCH --job-name=task12-empirical-holdout
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/evaluate_human_empirical_background_holdout.py \
  --candidate-model artifacts/task12_large_family_background/empirical_background_v1/model.joblib \
  --baseline-model artifacts/task12_large_family_background/pair_model_min20_q0.25_trim/model.joblib \
  --correction-reference artifacts/task10_human_calibration/replicate_0.seg \
  --holdouts \
    artifacts/task10_human_calibration/replicate_12.seg \
    artifacts/task10_human_calibration/replicate_13.seg \
    artifacts/task10_human_calibration/replicate_14.seg \
    artifacts/task10_human_calibration/replicate_15.seg \
  --out artifacts/task12_large_family_background/empirical_background_v1/holdout_summary.json
