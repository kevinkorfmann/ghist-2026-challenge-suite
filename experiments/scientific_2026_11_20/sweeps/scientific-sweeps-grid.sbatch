#!/usr/bin/env bash
#SBATCH --job-name=science-sweeps-grid
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-sweeps-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-sweeps-20261120/logs/%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/scientific-sweeps-20261120"
export LD_LIBRARY_PATH="$ROOT/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$PWD/src:$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 XDG_CACHE_HOME="$ROOT/.cache"
for OFFSET in $(seq 0 4 40); do
  INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))
  if ((INDEX < 20)); then
    COUNT=0
    SEED=$((202611300 + INDEX))
    OUTPUT="artifacts/grid/calibration_${INDEX}.json"
  else
    EVAL_INDEX=$((INDEX - 20))
    COUNT=$((EVAL_INDEX / 8))
    SEED=$((202611400 + EVAL_INDEX))
    OUTPUT="artifacts/grid/evaluation_${EVAL_INDEX}.json"
  fi
  "$ROOT/.venv/bin/python" -m scripts.validate_sweep_background_smoke --out "$OUTPUT" --coding-bed "$ROOT/data/reference/sweeps/havana104_chr7_exons.bed" --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --raisd "$ROOT/external/raisd/RAiSD" --sweeps "$COUNT" --seed "$SEED"
done
