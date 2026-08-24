#!/usr/bin/env bash
#SBATCH --job-name=task12-background-trim-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

quantiles=(0.25 0.50 0.75)
quantile="${quantiles[$SLURM_ARRAY_TASK_ID]}"
PYTHONPATH=src:scripts .venv/bin/python scripts/fit_human_family_relationships.py \
  --calibration \
    artifacts/task10_human_calibration/replicate_0.seg \
    artifacts/task10_human_calibration/replicate_1.seg \
    artifacts/task10_human_calibration/replicate_2.seg \
    artifacts/task10_human_calibration/replicate_3.seg \
    artifacts/task10_human_calibration/replicate_4.seg \
    artifacts/task10_human_calibration/replicate_5.seg \
  --validation artifacts/task10_human_calibration/replicate_6.seg \
  --evaluation artifacts/task10_human_calibration/replicate_7.seg \
  --observed artifacts/task12_large_family_background/hapibd/summary/pair_shape_features.csv \
  --background-reference artifacts/task10_large_family/hapibd/summary/pair_shape_features.csv \
  --background-quantile "$quantile" \
  --trim-background-shape \
  --minimum-cm 20 \
  --out-dir "artifacts/task12_large_family_background/pair_model_min20_q${quantile}_trim"
