#!/usr/bin/env bash
#SBATCH --job-name=ghist-task1-100mb
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-bottleneck-20261120/case-%A_%a.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/readiness-bottleneck-20261120"
export PYTHONPATH="$PWD/src:$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="$PWD/.mplconfig"
"$ROOT/.venv/bin/python" -m scripts.validate_bottleneck_readiness \
  --out "$PWD/results" --case "$SLURM_ARRAY_TASK_ID"
