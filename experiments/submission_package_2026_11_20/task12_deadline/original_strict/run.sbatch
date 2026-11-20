#!/bin/bash
#SBATCH --job-name=task12-final-fixed
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/final-task12-20261120/slurm-%j.log
source /etc/profile
set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/final-task12-20261120"
"$ROOT/.venv/bin/python" run.py --vcf "$ROOT/data/final-20261120/task12/relatedness4.final.vcf.gz" --binding "$ROOT/data/final-20261120/task12/input_binding.json" --out "$ROOT/scientific/final-task12-20261120"
