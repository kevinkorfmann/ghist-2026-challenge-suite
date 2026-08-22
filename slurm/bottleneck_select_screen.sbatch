#!/usr/bin/env bash
#SBATCH --job-name=bn-select-screen
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
SCREEN_JOB="${SCREEN_JOB:?set SCREEN_JOB to the completed screen array job}"
.venv/bin/python scripts/select_bottleneck_candidate.py \
  --root "checkpoints/bottleneck/screens/$SCREEN_JOB" \
  --out "artifacts/task1_bottleneck/screen_${SCREEN_JOB}_selection.json"
