#!/bin/bash
#SBATCH --job-name=ghist-lowcov-window-smoke
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/winsfs-smoke-%j.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
cd "$TASK_ROOT/source_v6"
export PYTHONPATH="$TASK_ROOT/source_v6/src:$TASK_ROOT/source_v6"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_winsfs.py \
 --case "$TASK_ROOT/smoke_v2_bottleneck" --output "$TASK_ROOT/window_smoke_bottleneck_v6" \
 --existing-inference "$TASK_ROOT/iid_gl2_bottleneck" --winsfs "$TASK_ROOT/external/winsfs/target/release/winsfs" --threads 4
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_winsfs.py \
 --case "$TASK_ROOT/smoke_v2_split_isolation" --output "$TASK_ROOT/window_smoke_split_isolation_v6" \
 --angsd "$TASK_ROOT/external/angsd/angsd" --winsfs "$TASK_ROOT/external/winsfs/target/release/winsfs" --threads 4
