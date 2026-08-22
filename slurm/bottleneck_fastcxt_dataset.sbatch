#!/usr/bin/env bash
#SBATCH --job-name=bn-fastcxt-data
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=24
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
DATASET_TAG="${DATASET_TAG:-fastcxt_v1}"
CELLS="${CELLS:-256}"
REPLICAS="${REPLICAS:-2}"
PAIRS="${PAIRS:-100}"

.venv/bin/python scripts/build_fastcxt_bottleneck_dataset.py \
  --out "data/simulations/bottleneck/$DATASET_TAG" \
  --cells "$CELLS" \
  --replicas "$REPLICAS" \
  --pairs "$PAIRS" \
  --workers 24
