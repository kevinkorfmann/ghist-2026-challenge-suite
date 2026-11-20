#!/bin/bash
#SBATCH --job-name=ghist-lowcov-matched-gl2
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/matched-gl2-%j.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
export PYTHONPATH="$TASK_ROOT/source_v5/src:$TASK_ROOT/source_v5"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_realsfs_from_saf_v5.py" \
 --case "$TASK_ROOT/smoke_v2_split_isolation" --output "$TASK_ROOT/matched_gl2_split_isolation_v5" \
 --saf-directory "$TASK_ROOT/window_smoke_split_isolation_v6/safs" \
 --realsfs "$TASK_ROOT/external/angsd/misc/realSFS" --threads 4
