#!/bin/bash
#SBATCH --job-name=ghist-lowcov-continuation
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/continue-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
cd "$TASK_ROOT/source_v4"
export PYTHONPATH="$TASK_ROOT/source_v4/src:$TASK_ROOT/source_v4"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
TASKS=(bottleneck split_isolation)
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_continue.py \
 --case "$TASK_ROOT/smoke_v2_${TASKS[$SLURM_ARRAY_TASK_ID]}" \
 --output "$TASK_ROOT/continuation_${TASKS[$SLURM_ARRAY_TASK_ID]}" --iterations 5000 --restarts 8
