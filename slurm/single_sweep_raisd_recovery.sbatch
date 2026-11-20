#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-raisd-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-63
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
INPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_recovery/task_${TASK_INDEX}"
OUTPUT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/selection_recovery_raisd/task_${TASK_INDEX}.json"
cd "$PROJECT_ROOT"

export LD_LIBRARY_PATH="$PROJECT_ROOT/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
PYTHONPATH=src .venv/bin/python scripts/run_single_sweep_raisd_recovery.py \
  --trees "$INPUT_ROOT.trees" \
  --truth "$INPUT_ROOT.json" \
  --raisd external/raisd/RAiSD \
  --out "$OUTPUT_ROOT"

