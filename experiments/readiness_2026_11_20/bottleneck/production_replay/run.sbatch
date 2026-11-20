#!/usr/bin/env bash
#SBATCH --job-name=ghist-production-replay
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/bottleneck-production-replay-20261120/logs/%A_%a.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
PROJECT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$PROJECT/reviews/bottleneck-production-replay-20261120"
export PYTHONPATH="$PWD/source_snapshot/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$PWD/.mplconfig"
"$PROJECT/.venv/bin/python" "$PWD/run_case.py" --case-index "$SLURM_ARRAY_TASK_ID"
