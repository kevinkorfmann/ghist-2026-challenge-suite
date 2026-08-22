#!/usr/bin/env bash
#SBATCH --job-name=bn-multi-eval
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/evaluate_multiscale_sfs_recovery.py \
  --input-glob 'data/simulations/bottleneck/multiscale_precision/multiscale_sfs_*.npz' \
  --observed artifacts/task1_bottleneck/moments_folded.json \
  --out artifacts/task1_bottleneck/multiscale_sfs_gate.json \
  --model-out artifacts/task1_bottleneck/multiscale_sfs_calibrator.joblib
