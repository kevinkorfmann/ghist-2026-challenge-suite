#!/bin/bash
#SBATCH --job-name=task10-comparator
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/comparator-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$D/complete_small_comparator.py"
