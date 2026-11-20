#!/bin/bash
#SBATCH --job-name=readiness-sweeps
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/readiness_sweeps/logs/%A_%a.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
CODE="$ROOT/reviews/readiness-sweeps-20261120"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME="$ROOT/.cache"
export LD_LIBRARY_PATH="$ROOT/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$CODE/src:$CODE"
cd "$CODE"
"$ROOT/.venv/bin/python" -m scripts.run_sweep_readiness_case \
  --design experiments/readiness_2026_11_20/sweeps/task5_design_v4.json \
  --case "$SLURM_ARRAY_TASK_ID" \
  --accessible "$ROOT/data/raw/shared/sweeps/testing.accessible.bed" \
  --slim "$ROOT/external/SLiM-v4.3/build/slim" \
  --raisd "$ROOT/external/raisd/RAiSD" \
  --out-dir "$ROOT/scientific/readiness_sweeps/task5_v4_memory_retry"
