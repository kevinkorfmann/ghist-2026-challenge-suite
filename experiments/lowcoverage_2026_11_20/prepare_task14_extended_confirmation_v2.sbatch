#!/bin/bash
#SBATCH --job-name=ghist-lowcov-freshprep
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --array=0-5%6
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/extend-v2-prep-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
export PYTHONPATH="$TASK_ROOT/source_v10/src:$TASK_ROOT/source_v10"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v10"
REGIME=$((SLURM_ARRAY_TASK_ID / 2))
REPLICATE=$((SLURM_ARRAY_TASK_ID % 2))
for_phase() {
 "$TASK_PYTHON" scripts/lowcoverage_clock_revision_guard.py --root "$TASK_ROOT" --freeze-sha "$FROZEN_SHA" \
 --regime "$REGIME" --replicate "$REPLICATE" --phase "$1"
}
for_phase prepare_initial
"$TASK_PYTHON" scripts/lowcoverage_prepare_confirmation.py --protocol "$TASK_ROOT/task14_extension_protocol_v2.json" \
 --task split_isolation --regime "$REGIME" --replicate "$REPLICATE" \
 --output "$TASK_ROOT/full_extended_confirmation_split_r${REGIME}_p${REPLICATE}_v1" \
 --bwa "$TASK_ROOT/external/bwa/bwa" --threads 4
for_phase prepare_final
