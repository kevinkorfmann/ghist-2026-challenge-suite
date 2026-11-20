#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --array=0-7%4
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/crossedfuture_%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/crossedfuture_%A_%a.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v5"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
FUTURE_SPLIT=$1
CANDIDATE_SHA=$2
if [[ "$FUTURE_SPLIT" == confirmation ]]; then
    "$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/run_crossed_future_v1.py" "$TASK_ROOT" run --split confirmation --index "$SLURM_ARRAY_TASK_ID" --candidate-sha256 "$CANDIDATE_SHA" --stress-summary-sha256 "$3"
elif [[ "$FUTURE_SPLIT" == stress ]]; then
    "$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/run_crossed_future_v1.py" "$TASK_ROOT" run --split stress --index "$SLURM_ARRAY_TASK_ID" --candidate-sha256 "$CANDIDATE_SHA"
else
    exit 2
fi
