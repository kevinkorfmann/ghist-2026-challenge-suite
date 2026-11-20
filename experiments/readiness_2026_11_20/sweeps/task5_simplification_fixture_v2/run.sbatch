#!/bin/bash
#SBATCH --job-name=sweep-simplify-fixture
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:05:00
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/sweep-simplification-fixture-v2-20261120"
"$ROOT/.venv/bin/python" run.py
