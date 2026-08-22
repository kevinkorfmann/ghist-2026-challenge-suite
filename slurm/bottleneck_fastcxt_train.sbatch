#!/usr/bin/env bash
#SBATCH --job-name=bn-fastcxt-train
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
DATASET_TAG="${DATASET_TAG:-fastcxt_v1}"
EPOCHS="${EPOCHS:-10}"

.venv/bin/fastcxt-train \
  --model small \
  --dataset-path "data/simulations/bottleneck/$DATASET_TAG" \
  --dataset-mode lazy \
  --gpus 0 \
  --epochs "$EPOCHS" \
  --batch-size 128 \
  --grad-accum 4 \
  --workers 6 \
  --no-compile \
  --csv-log \
  --save-every-epoch \
  --log-dir "checkpoints/fastcxt/bottleneck/${DATASET_TAG}_${SLURM_JOB_ID}"
