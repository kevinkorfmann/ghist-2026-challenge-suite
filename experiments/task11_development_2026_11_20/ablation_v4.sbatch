#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v4"
export PYTHONPATH="$PWD:$PWD/src:$PWD/scripts"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT_ROOT/.venv/bin/python" "$TASK_ROOT/phase_development_ablation.py" "$TASK_ROOT"
