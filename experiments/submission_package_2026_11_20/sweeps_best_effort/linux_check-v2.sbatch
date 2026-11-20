#!/bin/bash
#SBATCH --job-name=sweep-final-check
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:03:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120/linux_checks_v2.log
source /etc/profile >/dev/null 2>&1
module load python/3.11.11
set -euo pipefail
cd /vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120
export PYTHONPATH="$PWD/src:$PWD/scripts:/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/2026-11-20/.test-deps"
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python -m pytest --rootdir="$PWD" --confcutdir="$PWD" -q experiments/submission_package_2026_11_20/sweeps_best_effort/test_runner.py
