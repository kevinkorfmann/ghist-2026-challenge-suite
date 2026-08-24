#!/usr/bin/env bash
#SBATCH --job-name=task12-background-robust
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/fit_human_background_robust.py \
  --calibration \
    artifacts/task10_human_calibration/replicate_0.seg \
    artifacts/task10_human_calibration/replicate_1.seg \
    artifacts/task10_human_calibration/replicate_2.seg \
    artifacts/task10_human_calibration/replicate_3.seg \
    artifacts/task10_human_calibration/replicate_4.seg \
    artifacts/task10_human_calibration/replicate_5.seg \
  --selection artifacts/task10_human_calibration/replicate_6.seg \
  --development-evaluation artifacts/task10_human_calibration/replicate_7.seg \
  --observed artifacts/task12_large_family_background/hapibd/summary/pair_shape_features.csv \
  --minimum-cm 20 \
  --augmentations 4 \
  --out-dir artifacts/task12_large_family_background/background_robust_v1
