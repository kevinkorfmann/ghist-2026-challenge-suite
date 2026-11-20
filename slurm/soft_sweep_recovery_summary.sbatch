#!/usr/bin/env bash
#SBATCH --job-name=soft-sweep-summary
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/summarize_soft_sweep_recovery.py \
  --input-dir artifacts/task6_soft_sweep/recovery_feasible \
  --expected 48 \
  --out artifacts/task6_soft_sweep/recovery_summary.json
