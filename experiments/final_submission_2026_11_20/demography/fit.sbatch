#!/usr/bin/env bash
#SBATCH --job-name=final-demography-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem=5632M
#SBATCH --time=00:20:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/final-demography-20261120/logs/fit-%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/final-demography-20261120/logs/fit-%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
cd /vast/projects/smathi/cohort/korfmann/ghist2026/reviews/final-demography-20261120
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python run.py fit "$SLURM_ARRAY_TASK_ID"
