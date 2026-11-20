#!/usr/bin/env bash
#SBATCH --job-name=science-sweeps-blind
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
export PYTHONPATH="$PWD/src:$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 XDG_CACHE_HOME="$ROOT/.cache"
for OFFSET in $(seq 0 4 188); do
  INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))
  INPUT="$ROOT/artifacts/task5_single_sweep/selection_exact_recovery/task_${INDEX}"
  "$ROOT/.venv/bin/python" -m scripts.recompute_sweep_recovery --record "$INPUT.json" --trees "$INPUT.trees" --scan "$INPUT.csv" --out "artifacts/blind/task_${INDEX}.json"
done
