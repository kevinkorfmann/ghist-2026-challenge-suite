#!/usr/bin/env bash
#SBATCH --job-name=human-ibd-shape
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  root="artifacts/task10_large_family/hapibd"
else
  root="artifacts/task12_large_family_background/hapibd"
fi

PYTHONPATH=src:scripts .venv/bin/python scripts/summarize_hapibd_shape.py \
  --ibd "$root/segments.ibd.gz" \
  --base "$root/summary/pair_features.csv" \
  --out "$root/summary/pair_shape_features.csv"
