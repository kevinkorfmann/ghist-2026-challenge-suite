#!/bin/bash
#SBATCH --job-name=ghist-lowcov-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/audit-spectra-%j.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
export PYTHONPATH="$TASK_ROOT/source_v5/src:$TASK_ROOT/source_v5"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_audit_spectra_v5.py" \
 --case "$TASK_ROOT/smoke_v2_bottleneck" --output "$TASK_ROOT/truth_control_bottleneck_v5"
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_audit_spectra_v5.py" \
 --case "$TASK_ROOT/smoke_v2_bottleneck" --output "$TASK_ROOT/gl1_refit_bottleneck_v5" \
 --sfs "$TASK_ROOT/smoke_v2_bottleneck/inference/folded.sfs"
