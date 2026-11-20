#!/bin/bash
#SBATCH --job-name=task10-dev-diagnose
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/diagnostic-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
S="$A/experiments/task10-development-20261120/diagnostic-source"
export PYTHONPATH="$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$A/.venv/bin/python" "$S/scripts/develop_task10_detected_ibd.py" diagnose --cohort "$A/artifacts/task10-complete-20261120/seed_2026112010" --asset-root "$A" --out-dir "$A/experiments/task10-development-20261120/diagnostic"
