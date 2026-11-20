#!/bin/bash
#SBATCH --job-name=task10-site-thin
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/site-thin-%A_%a.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/structural-qualification-source-v3"
H="$D/site-thinning-source-v1"
export PYTHONPATH="$Q:$S/experiments/task10_development_2026_11_20:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$H/run_site_thinning_diagnostic.py" --source "$H" --root "$D" --index "$SLURM_ARRAY_TASK_ID"
