#!/usr/bin/env bash
#SBATCH --job-name=human-family-calibration
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
PED_ROOT="$PROJECT_ROOT/external/ped-sim"
MAP_ROOT="$PROJECT_ROOT/external/ped-sim-map"
OUT_ROOT="$PROJECT_ROOT/artifacts/task10_human_calibration"
mkdir -p "$OUT_ROOT"

seed=$((101226000 + SLURM_ARRAY_TASK_ID))
"$PED_ROOT/ped-sim" \
  -d "$PROJECT_ROOT/configs/task10_human_pedigrees.def" \
  -m "$MAP_ROOT/refined_mf.simmap" \
  -o "$OUT_ROOT/replicate_${SLURM_ARRAY_TASK_ID}" \
  --intf "$PED_ROOT/interfere/nu_p_campbell.tsv" \
  --seed "$seed"
