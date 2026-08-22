#!/usr/bin/env bash
#SBATCH --job-name=bn-screen
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --array=0-4
#SBATCH --cpus-per-task=6
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
CANDIDATES=(baseline wider deeper_context extended_local compact)
CANDIDATE="${CANDIDATES[$SLURM_ARRAY_TASK_ID]}"
.venv/bin/python scripts/run_bottleneck_experiment.py \
  --registry configs/bottleneck_experiments.yaml \
  --candidate "$CANDIDATE" \
  --dataset data/simulations/bottleneck/search.npz \
  --out-root "checkpoints/bottleneck/screens/${SLURM_ARRAY_JOB_ID}"
