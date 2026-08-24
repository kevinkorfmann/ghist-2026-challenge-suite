#!/usr/bin/env bash
#SBATCH --job-name=task11-ibd-holdout
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=8-15
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
PED_ROOT="$PROJECT_ROOT/external/ped-sim"
MAP_ROOT="$PROJECT_ROOT/external/dog-recombination"
OUT_ROOT="$PROJECT_ROOT/artifacts/task11_multifamily_inbreeding/calibration"
mkdir -p "$OUT_ROOT"

seed=$((91126000 + SLURM_ARRAY_TASK_ID))
"$PED_ROOT/ped-sim" \
  -d "$PROJECT_ROOT/configs/task11_pedigrees.def" \
  -m "$MAP_ROOT/campbell2016_average.simmap" \
  -o "$OUT_ROOT/replicate_${SLURM_ARRAY_TASK_ID}" \
  --pois \
  --seed "$seed"
