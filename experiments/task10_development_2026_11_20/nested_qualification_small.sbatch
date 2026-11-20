#!/bin/bash
#SBATCH --job-name=task10-nested-small
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/nested-small-%A_%a.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/nested-qualification-source-v1"
P="$S/experiments/task10_development_2026_11_20"
export PYTHONPATH="$Q:$P:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/usr/bin/time -v "$A/.venv/bin/python" "$Q/run_nested_qualification.py" --frozen-source "$S" --qualification-source "$Q" --out-root "$D/nested_qualification_v1" --index "$SLURM_ARRAY_TASK_ID"
