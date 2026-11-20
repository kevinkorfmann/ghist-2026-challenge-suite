#!/usr/bin/env bash
#SBATCH --job-name=relatedness-missingness
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=0-3
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/relatedness-missingness-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/relatedness-missingness-20261120/logs/%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
TASK_ROOT="$PROJECT_ROOT/experiments/relatedness-missingness-20261120"
cd "$TASK_ROOT/source"
export PYTHONPATH="$PWD:$PWD/src:$PWD/scripts"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
seeds=(2026112711 2026112811)
rates=(0.005 0.02)
seed="${seeds[$((SLURM_ARRAY_TASK_ID % 2))]}"
rate="${rates[$((SLURM_ARRAY_TASK_ID / 2))]}"
"$PROJECT_ROOT/.venv/bin/python" scripts/validate_relatedness_missingness.py \
  --baseline "$PROJECT_ROOT/experiments/task11-complete-20261120/results/seed_$seed" \
  --rate "$rate" \
  --out "$TASK_ROOT/results/seed_${seed}_rate_${rate}"
