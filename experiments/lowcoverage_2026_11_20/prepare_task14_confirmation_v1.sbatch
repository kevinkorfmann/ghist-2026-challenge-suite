#!/bin/bash
#SBATCH --job-name=ghist-lowcov-confirm-prep
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=0-5%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/confirm-prep-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
REGIME=$((SLURM_ARRAY_TASK_ID / 2))
REPLICATE=$((SLURM_ARRAY_TASK_ID % 2))
STEM=split_r${REGIME}_p${REPLICATE}_v1
export PYTHONPATH="$TASK_ROOT/source_v7/src:$TASK_ROOT/source_v7"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v7"
guard() {
 "$TASK_PYTHON" "$TASK_ROOT/lowcoverage_confirm_guard_v1.py" --root "$TASK_ROOT" \
 --freeze-sha "$FROZEN_SHA" --regime "$REGIME" --replicate "$REPLICATE" --phase "$1"
}
guard prepare_initial
"$TASK_PYTHON" scripts/lowcoverage_prepare_confirmation.py --protocol "$TASK_ROOT/validation_acceptance_v1.json" \
 --task split_isolation --regime "$REGIME" --replicate "$REPLICATE" --output "$TASK_ROOT/full_confirmation_$STEM" \
 --bwa "$TASK_ROOT/external/bwa/bwa" --threads 4
guard prepare_final
