#!/bin/bash
#SBATCH --job-name=ghist-task2-regression
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=logs/task2-regression-%j.out
#SBATCH --error=logs/task2-regression-%j.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/2026-11-20"
export PYTHONPATH="$PWD/src:$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$ROOT/.venv/bin/python" artifacts/reliability/scan_task2.py
