#!/usr/bin/env bash
#SBATCH --job-name=bn-merge-search
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/merge_simulations.py \
  --input 'data/simulations/bottleneck/part_*.npz' \
  --output data/simulations/bottleneck/search.npz
