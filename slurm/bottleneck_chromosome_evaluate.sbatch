#!/usr/bin/env bash
#SBATCH --job-name=bn-chrom-eval
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
PROMOTION_JOB="${PROMOTION_JOB:?set PROMOTION_JOB to the completed promotion job}"
.venv/bin/python scripts/evaluate_bottleneck_chromosomes.py \
  --dataset data/simulations/bottleneck/chromosome_validation.npz \
  --checkpoint-dir "checkpoints/bottleneck/promotions/fresh_${PROMOTION_JOB}/deeper_context" \
  --out "artifacts/task1_bottleneck/chromosome_validation_${PROMOTION_JOB}.json"
