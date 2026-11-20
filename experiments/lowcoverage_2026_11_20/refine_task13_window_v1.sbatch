#!/bin/bash
#SBATCH --job-name=ghist-lowcov-refine
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --array=0-2%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/refine-task13-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
LABELS=(seed1 seed2 continued)
LABEL=${LABELS[$SLURM_ARRAY_TASK_ID]}
export PYTHONPATH="$TASK_ROOT/source_v7/src:$TASK_ROOT/source_v7"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v7"
"$TASK_PYTHON" "$TASK_ROOT/lowcoverage_refine_window_v1.py" --root "$TASK_ROOT" --label "$LABEL" \
 --output "$TASK_ROOT/full_refinement_bottleneck_${LABEL}_v1"
