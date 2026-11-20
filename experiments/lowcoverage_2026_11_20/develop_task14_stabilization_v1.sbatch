#!/bin/bash
#SBATCH --job-name=ghist-lowcov-stabilize
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=0-1%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/stabilize-development-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
export PYTHONPATH="$TASK_ROOT/source_v8/src:$TASK_ROOT/source_v8"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v8"
"$TASK_PYTHON" scripts/lowcoverage_stabilize.py --root "$TASK_ROOT" --freeze-sha "$FROZEN_SHA" \
 --mode development --regime 0 --replicate "$SLURM_ARRAY_TASK_ID"
