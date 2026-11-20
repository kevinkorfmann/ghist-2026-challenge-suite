#!/bin/bash
#SBATCH --job-name=ghist-lowcov-full-prep
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --array=0-1
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/full-prep-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
cd "$TASK_ROOT/source_v3"
export PYTHONPATH="$TASK_ROOT/source_v3/src:$TASK_ROOT/source_v3"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
TASKS=(bottleneck split_isolation)
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_prepare_full.py \
 --task "${TASKS[$SLURM_ARRAY_TASK_ID]}" --seed-offset 10000 \
 --output "$TASK_ROOT/full_${TASKS[$SLURM_ARRAY_TASK_ID]}" \
 --bwa "$TASK_ROOT/external/bwa/bwa" --threads 4
