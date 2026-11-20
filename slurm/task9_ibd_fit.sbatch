#!/usr/bin/env bash
#SBATCH --job-name=task9-ibd-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/fit_task9_ibd_classifier.py \
  --calibration \
    artifacts/task9_ibd_segments/calibration/replicate_0.seg \
    artifacts/task9_ibd_segments/calibration/replicate_1.seg \
    artifacts/task9_ibd_segments/calibration/replicate_2.seg \
    artifacts/task9_ibd_segments/calibration/replicate_3.seg \
    artifacts/task9_ibd_segments/calibration/replicate_4.seg \
  --validation artifacts/task9_ibd_segments/calibration/replicate_5.seg \
  --evaluation artifacts/task9_ibd_segments/calibration/replicate_6.seg \
  --observed data/raw/relatedness/ibd_segments/testing/relatedness1.testing.csv \
  --sex-map external/ped-sim-map/refined_mf.simmap \
  --out-dir artifacts/task9_ibd_segments/model_v2
