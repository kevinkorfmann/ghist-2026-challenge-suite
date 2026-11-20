#!/bin/bash
#SBATCH --job-name=sweep-final
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
source /etc/profile >/dev/null 2>&1
module load python/3.11.11
set -euo pipefail
SWEEP_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120
export PYTHONPATH="$SWEEP_ROOT/src:$SWEEP_ROOT/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
exec /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/run_final_sweep.py" --source-root "$SWEEP_ROOT" --source-freeze "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/source_freeze.json" --raisd /vast/projects/smathi/cohort/korfmann/ghist2026/external/raisd/RAiSD "$@"
