#!/usr/bin/env bash
#SBATCH --job-name=bn-promotion
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
RUN_ROOT="checkpoints/bottleneck/promotions/fresh_${SLURM_JOB_ID}"

.venv/bin/python scripts/run_bottleneck_experiment.py \
  --registry configs/bottleneck_promotion.yaml \
  --candidate deeper_context \
  --dataset data/simulations/bottleneck/promotion.npz \
  --out-root "$RUN_ROOT"

.venv/bin/python -m ghist_bottleneck.cli infer-blocks \
  --observed data/processed/task1_bottleneck/observed.blocks.npz \
  --checkpoint-dir "$RUN_ROOT/deeper_context" \
  --out "artifacts/task1_bottleneck/promotion_${SLURM_JOB_ID}_summary.yaml"

.venv/bin/python -m ghist_bottleneck.cli infer-blocks \
  --observed data/processed/task1_bottleneck/observed.fastcxt.blocks.npz \
  --checkpoint-dir "$RUN_ROOT/deeper_context" \
  --out "artifacts/task1_bottleneck/promotion_${SLURM_JOB_ID}_fastcxt.yaml"
