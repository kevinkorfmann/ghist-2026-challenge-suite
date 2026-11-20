#!/bin/bash
#SBATCH --job-name=task10-tree-all-development
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --array=0-2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/tree-all-%A_%a.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/tree-uncertainty-source"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
case "$SLURM_ARRAY_TASK_ID" in
0) C="$D/diagnostic"; T="$A/artifacts/task10-complete-20261120/seed_2026112010/truth.json";;
1) C="$D/development_cohort_0/representation_comparison"; T="$D/development_cohort_0/truth.json";;
2) C="$D/development_cohort_1/representation_comparison"; T="$D/development_cohort_1/truth.json";;
esac
"$A/.venv/bin/python" "$S/scripts/develop_task10_detected_ibd.py" tree --cohort "$C" --truth "$T" --asset-root "$A" --out-dir "$D/tree_analysis_$SLURM_ARRAY_TASK_ID"
