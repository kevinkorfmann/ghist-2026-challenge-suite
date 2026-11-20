#!/bin/bash
#SBATCH --job-name=task7-software
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/readiness_sweeps/logs/task7-software-%j.log
source /etc/profile
set -euo pipefail
module load python/3.11.11
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
SNAPSHOT="$ROOT/reviews/task7-software-20261120"
export PYTHONPATH="$SNAPSHOT/src:$SNAPSHOT"
export LD_LIBRARY_PATH="$ROOT/external/raisd/gsl/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$SNAPSHOT"
"$ROOT/.venv/bin/python" run.py "$ROOT/scientific/readiness_sweeps/task7_software_smoke" "$ROOT/external/raisd/RAiSD"
