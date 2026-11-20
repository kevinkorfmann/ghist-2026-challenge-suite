#!/usr/bin/env bash
#SBATCH --job-name=splitmig-block-summary
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

PYTHONPATH=src .venv/bin/python scripts/summarize_split_migration_bgs_bootstrap.py \
  --task-dir artifacts/task4_split_migration_bgs/bootstrap_1mb \
  --block-summary artifacts/task4_split_migration_bgs/blocks_1mb_summary.json \
  --out artifacts/task4_split_migration_bgs/bootstrap_1mb_summary.json
