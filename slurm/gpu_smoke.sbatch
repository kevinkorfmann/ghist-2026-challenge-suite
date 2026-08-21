#!/usr/bin/env bash
#SBATCH --job-name=gpu-smoke
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:45gb:1
#SBATCH --cpus-per-task=6
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/gpu_smoke.py

