#!/usr/bin/env bash
#SBATCH --job-name=ghist-science-tests
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-integration-20261120/pytest-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/scientific-integration-20261120"
export PYTHONPATH="$ROOT/reviews/2026-11-20/.test-deps:$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/.mplconfig"
"$ROOT/.venv/bin/python" -m pytest -q
