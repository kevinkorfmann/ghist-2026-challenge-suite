#!/usr/bin/env bash
#SBATCH --job-name=task12-corrected-holdout
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/evaluate_human_background_corrected_holdout.py \
  --candidate-model artifacts/task12_large_family_background/background_corrected_min5_v2/model.joblib \
  --baseline-model artifacts/task12_large_family_background/pair_model_min20_q0.25_trim/model.joblib \
  --holdouts \
    artifacts/task10_human_calibration/replicate_8.seg \
    artifacts/task10_human_calibration/replicate_9.seg \
    artifacts/task10_human_calibration/replicate_10.seg \
    artifacts/task10_human_calibration/replicate_11.seg \
  --out artifacts/task12_large_family_background/background_corrected_min5_v2/holdout_summary.json
