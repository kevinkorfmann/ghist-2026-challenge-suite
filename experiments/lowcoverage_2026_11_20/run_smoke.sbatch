#!/bin/bash
#SBATCH --job-name=ghist-lowcov-smoke
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --array=0-1
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/smoke-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
cd "$TASK_ROOT/source_v2"
export PYTHONPATH="$TASK_ROOT/source_v2/src:$TASK_ROOT/source_v2_v2"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
TASKS=(bottleneck split_isolation)
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_validate.py \
 --task "${TASKS[$SLURM_ARRAY_TASK_ID]}" --length 1000000 \
 --output "$TASK_ROOT/smoke_v2_${TASKS[$SLURM_ARRAY_TASK_ID]}" \
 --angsd "$TASK_ROOT/external/angsd/angsd" --realsfs "$TASK_ROOT/external/angsd/misc/realSFS" \
 --threads 4 --restarts 8
