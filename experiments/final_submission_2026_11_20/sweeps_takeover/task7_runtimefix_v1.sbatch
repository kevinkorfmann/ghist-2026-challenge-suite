#!/bin/bash
#SBATCH --job-name=sweep-final7-runtimefix
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120/task7-runtimefix-%j.log
source /etc/profile >/dev/null 2>&1
module load python/3.11.11
set -euo pipefail
SWEEP_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120
SWEEP_DATA=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120
export LD_LIBRARY_PATH=/vast/projects/smathi/cohort/korfmann/ghist2026/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$SWEEP_ROOT/src:$SWEEP_ROOT/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
exec /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/run_final_sweep.py" --source-root "$SWEEP_ROOT" --source-freeze "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/source_freeze.json" --raisd /vast/projects/smathi/cohort/korfmann/ghist2026/external/raisd/RAiSD --task 7 --vcf "$SWEEP_DATA/inputs/task7/onesweep_bgs.final.vcf.gz" --vcf-sha256 ef4f6fb9a795df476f6553e5890db1b3094fc026f288f11378b6c05c1ef67baa --accessible "$SWEEP_DATA/inputs/shared/9.accessible.bed" --mask-sha256 9f5235a687683eec6a7f19b0cf1334f545bda23497ff8e715658b1086778ec33 --out-dir "$SWEEP_DATA/task7_runtimefix_v1"
