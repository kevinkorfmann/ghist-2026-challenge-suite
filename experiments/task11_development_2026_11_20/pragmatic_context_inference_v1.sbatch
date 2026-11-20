#!/usr/bin/env bash
# Root-owned inference template; preparing this file does not dispatch it.
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v5"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/prepare_pragmatic_context_v1.py" "$TASK_ROOT" \
  --vcf "$1" --out "$2" --dataset final \
  --candidate-sha256 6ccfda5f8133b5dee5265f35ac3a328fc4c35cfa88af076e11b9873af781ea15 \
  --runner-sha256 7e2bdb9ffcd8f42802e2470a93f2c325b6857d57a963854b74b37a41ef1f04ce
