#!/usr/bin/env bash
#SBATCH --job-name=ghist-case08-numeric
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/case08-numerical-20261120/slurm-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
PROJECT=/vast/projects/smathi/cohort/korfmann/ghist2026
DIAGNOSTIC="$PROJECT/reviews/case08-numerical-20261120"
cd "$DIAGNOSTIC"
export PYTHONPATH="$PROJECT/reviews/readiness-bottleneck-20261120/src:$PROJECT/reviews/readiness-bottleneck-20261120"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$DIAGNOSTIC/.mplconfig"
"$PROJECT/.venv/bin/python" "$DIAGNOSTIC/diagnose.py" \
  --case "$DIAGNOSTIC/input" --out "$DIAGNOSTIC/results"
