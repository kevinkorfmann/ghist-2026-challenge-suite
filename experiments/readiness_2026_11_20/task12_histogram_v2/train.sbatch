#!/bin/bash
#SBATCH --job-name=task12-histogram-train
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task12-histogram-v2-20261120/train-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task12-histogram-v2-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TASK_ROOT/.venv/bin/python" scripts/train_task12_histogram_candidate.py --root "$TASK_ROOT" \
  --protocol experiments/readiness_2026_11_20/task12_histogram_v2/frozen_protocol.json \
  --out-dir artifacts/task12_histogram_v2_20261120
