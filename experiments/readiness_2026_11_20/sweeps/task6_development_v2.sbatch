#!/bin/bash
#SBATCH --job-name=readiness-soft-dev
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-47%4
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/readiness_sweeps/logs/softdevv2_%A_%a.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
CODE="$ROOT/reviews/readiness-sweeps-development-v2-20261120"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$CODE/src:$CODE"
cd "$CODE"
"$ROOT/.venv/bin/python" -m scripts.develop_soft_sweep \
 --trees "$ROOT/artifacts/task6_soft_sweep/recovery_feasible/task_${SLURM_ARRAY_TASK_ID}.trees" \
 --record "$ROOT/artifacts/task6_soft_sweep/recovery_feasible/task_${SLURM_ARRAY_TASK_ID}.json" \
 --out "$ROOT/scientific/readiness_sweeps/task6_development_v2/task_${SLURM_ARRAY_TASK_ID}.json"
