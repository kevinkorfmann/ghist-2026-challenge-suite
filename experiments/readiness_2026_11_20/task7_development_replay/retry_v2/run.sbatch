#!/usr/bin/env bash
#SBATCH --job-name=task7-dev-replay-v2
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/task7-development-replay-v2-20261120/slurm-%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/task7-development-replay-v2-20261120/slurm-%j.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
REPLAY_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/task7-development-replay-v2-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$REPLAY_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export LD_LIBRARY_PATH="$PROJECT_ROOT/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$REPLAY_ROOT/source_snapshot/src:$REPLAY_ROOT/source_snapshot"
"$PROJECT_ROOT/.venv/bin/python" run.py
