#!/usr/bin/env bash
#SBATCH --job-name=task10-spectrum-holdout
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/evaluate_human_clean_spectrum_holdout.py \
  --candidate-model artifacts/task10_large_family/clean_spectrum_v1/model.joblib \
  --baseline-model artifacts/task10_large_family/pair_model_min5/model.joblib \
  --holdouts \
    artifacts/task10_human_calibration/replicate_16.seg \
    artifacts/task10_human_calibration/replicate_17.seg \
    artifacts/task10_human_calibration/replicate_18.seg \
    artifacts/task10_human_calibration/replicate_19.seg \
  --out artifacts/task10_large_family/clean_spectrum_v1/holdout_summary.json
