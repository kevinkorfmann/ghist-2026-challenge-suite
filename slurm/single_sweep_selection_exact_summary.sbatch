#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-exact-summary
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/summarize_single_sweep_exact_recovery.py \
  --input-dir artifacts/task5_single_sweep/selection_exact_recovery \
  --expected 192 \
  --out artifacts/task5_single_sweep/selection_exact_recovery_summary.json
