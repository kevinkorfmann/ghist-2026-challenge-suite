#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-extract
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
WINDOW_INDEX="${WINDOW_INDEX:-500}"
RUN_ROOT="$PROJECT_ROOT/artifacts/task2_secondary_contact/singer/window_${WINDOW_INDEX}"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/extract_secondary_contact_arg_features.py \
  --trees "$RUN_ROOT/topology_*.trees" \
  --out "$RUN_ROOT/coalescence_hazards.npz"
