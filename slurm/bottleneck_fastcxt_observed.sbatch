#!/usr/bin/env bash
#SBATCH --job-name=bn-fastcxt-observed
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
TRAIN_JOB="${TRAIN_JOB:?set TRAIN_JOB to the completed FastCXT training job}"
DATASET_TAG="${DATASET_TAG:-fastcxt_v1}"

.venv/bin/python scripts/run_fastcxt_observed.py \
  --matrix-dir data/processed/task1_bottleneck/fastcxt_input \
  --checkpoint-root "checkpoints/fastcxt/bottleneck/${DATASET_TAG}_${TRAIN_JOB}" \
  --out artifacts/task1_bottleneck/fastcxt_observed.npz \
  --pairs 128

.venv/bin/python scripts/inject_fastcxt_hazards.py \
  --collection data/processed/task1_bottleneck/observed.blocks.npz \
  --translation artifacts/task1_bottleneck/fastcxt_observed.npz \
  --out data/processed/task1_bottleneck/observed.fastcxt.blocks.npz
