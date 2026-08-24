#!/usr/bin/env bash
#SBATCH --job-name=human-family-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=03:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

common=(
  --calibration
    artifacts/task10_human_calibration/replicate_0.seg
    artifacts/task10_human_calibration/replicate_1.seg
    artifacts/task10_human_calibration/replicate_2.seg
    artifacts/task10_human_calibration/replicate_3.seg
    artifacts/task10_human_calibration/replicate_4.seg
    artifacts/task10_human_calibration/replicate_5.seg
  --validation artifacts/task10_human_calibration/replicate_6.seg
  --evaluation artifacts/task10_human_calibration/replicate_7.seg
)

case "$SLURM_ARRAY_TASK_ID" in
  0|1|2)
    minimums=(5 10 20)
    minimum="${minimums[$SLURM_ARRAY_TASK_ID]}"
    observed="artifacts/task10_large_family/hapibd/summary/pair_shape_features.csv"
    out="artifacts/task10_large_family/pair_model_min${minimum}"
    extra=()
    ;;
  3|4|5)
    quantiles=(0.25 0.50 0.75)
    quantile="${quantiles[$((SLURM_ARRAY_TASK_ID - 3))]}"
    minimum=20
    observed="artifacts/task12_large_family_background/hapibd/summary/pair_shape_features.csv"
    out="artifacts/task12_large_family_background/pair_model_min20_q${quantile}"
    extra=(
      --background-reference artifacts/task10_large_family/hapibd/summary/pair_shape_features.csv
      --background-quantile "$quantile"
    )
    ;;
  6)
    minimum=40
    observed="artifacts/task12_large_family_background/hapibd/summary/pair_shape_features.csv"
    out="artifacts/task12_large_family_background/pair_model_min40_q0.50"
    extra=(
      --background-reference artifacts/task10_large_family/hapibd/summary/pair_shape_features.csv
      --background-quantile 0.50
    )
    ;;
esac

PYTHONPATH=src:scripts .venv/bin/python scripts/fit_human_family_relationships.py \
  "${common[@]}" \
  --observed "$observed" \
  --minimum-cm "$minimum" \
  --out-dir "$out" \
  "${extra[@]}"
