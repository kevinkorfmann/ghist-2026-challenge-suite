#!/bin/bash
#SBATCH --job-name=sweep-interval100-case5
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
CODE="$ROOT/reviews/readiness-sweeps-interval100-20261120"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME="$ROOT/.cache"
export LD_LIBRARY_PATH="$ROOT/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$CODE/snapshot/src:$CODE/snapshot"
cd "$CODE"
"$ROOT/.venv/bin/python" adapter.py
