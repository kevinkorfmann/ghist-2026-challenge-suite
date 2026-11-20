#!/bin/bash
#SBATCH --job-name=ghist-lowcov-gl2
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/gl2-smoke-%j.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
cd "$TASK_ROOT/source_v5"
export PYTHONPATH="$TASK_ROOT/source_v5/src:$TASK_ROOT/source_v5"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_reinfer.py \
 --case "$TASK_ROOT/smoke_v2_bottleneck" --output "$TASK_ROOT/iid_gl2_bottleneck" --gl-model 2 \
 --angsd "$TASK_ROOT/external/angsd/angsd" --realsfs "$TASK_ROOT/external/angsd/misc/realSFS" --threads 4
