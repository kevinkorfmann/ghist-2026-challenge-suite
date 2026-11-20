#!/usr/bin/env bash
#SBATCH --job-name=human-background-validation
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/validate_human_background_correction.py \
  --model artifacts/task12_large_family_background/pair_model_min20_q0.50/model.joblib \
  --evaluation artifacts/task10_human_calibration/replicate_7.seg \
  --reference artifacts/task10_human_calibration/replicate_0.seg \
  --replicates 12 \
  --out artifacts/task12_large_family_background/background_validation/summary.json
