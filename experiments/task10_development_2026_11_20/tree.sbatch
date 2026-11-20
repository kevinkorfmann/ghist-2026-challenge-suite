#!/bin/bash
#SBATCH --job-name=task10-tree-development
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/tree-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/tree-source"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$S"
"$A/.venv/bin/python" -m pytest tests/test_task10_detected_ibd.py -q
"$A/.venv/bin/python" "$S/scripts/develop_task10_detected_ibd.py" tree --cohort "$D/diagnostic" --truth "$A/artifacts/task10-complete-20261120/seed_2026112010/truth.json" --asset-root "$A" --out-dir "$D/tree_diagnostic"
