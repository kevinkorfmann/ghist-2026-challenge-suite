#!/usr/bin/env bash
#SBATCH --job-name=task12-deconvolution
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/fit_human_background_deconvolution.py \
  --model artifacts/task12_large_family_background/pair_model_min20_q0.25_trim/model.joblib \
  --selection artifacts/task10_human_calibration/replicate_6.seg \
  --development-evaluation artifacts/task10_human_calibration/replicate_7.seg \
  --observed artifacts/task12_large_family_background/hapibd/summary/pair_shape_features.csv \
  --out-dir artifacts/task12_large_family_background/background_deconvolution_v1
