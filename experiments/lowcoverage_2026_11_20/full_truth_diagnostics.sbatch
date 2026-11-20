#!/bin/bash
#SBATCH --job-name=ghist-lowcov-fulltruth
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/full-truth-%j.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
export PYTHONPATH="$TASK_ROOT/source_v5/src:$TASK_ROOT/source_v5"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
for TASK in bottleneck split_isolation; do
 /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_audit_spectra_v5.py" \
  --case "$TASK_ROOT/full_$TASK" --output "$TASK_ROOT/full_truth_${TASK}_v5"
done
for CASE in smoke_v2_bottleneck full_bottleneck; do
 /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_profile_truth_v5.py" \
  --case "$TASK_ROOT/$CASE" --output "$TASK_ROOT/profile_${CASE}_v5"
done
