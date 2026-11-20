#!/usr/bin/env bash
#SBATCH --job-name=science-sweeps
#SBATCH --partition=genoa-std-mem
#SBATCH --array=0-2
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
"$ROOT/.venv/bin/python" -m scripts.validate_sweep_background_smoke --out "artifacts/background_${SLURM_ARRAY_TASK_ID}.json" --coding-bed "$ROOT/data/reference/sweeps/havana104_chr7_exons.bed" --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --raisd "$ROOT/external/raisd/RAiSD" --sweeps "$SLURM_ARRAY_TASK_ID" --seed "$((202611200 + SLURM_ARRAY_TASK_ID))"
