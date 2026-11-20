#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/pragmatic_corrected_freeze_smoke_%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/pragmatic_corrected_freeze_smoke_%j.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v5"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/pragmatic_corrected_bindings_v1.py" "$TASK_ROOT" freeze --role base
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/pragmatic_corrected_bindings_v1.py" "$TASK_ROOT" freeze --role context
CANDIDATE_SHA=$(sha256sum "$TASK_ROOT/pragmatic_corrected_context_candidate_v1.json" | cut -d' ' -f1)
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/smoke_pragmatic_corrected_v1.py" "$TASK_ROOT" "$CANDIDATE_SHA" context
