#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-aggregate
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/aggregate_secondary_contact_arg_tiles.py \
  --root artifacts/task2_secondary_contact/singer_tiles_100kb \
  --out artifacts/task2_secondary_contact/tiled_arg_features.npz \
  --tiles 1000 --tiles-per-block 20
