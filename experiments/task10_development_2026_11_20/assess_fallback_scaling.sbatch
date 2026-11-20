#!/bin/bash
#SBATCH --job-name=task10-fresh-assess
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/fallback-assess-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/structural-qualification-source-v3"
E="$D/fallback-assessment-source-v1"
export PYTHONPATH="$E:$Q:$S/experiments/task10_development_2026_11_20:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$E/collect_structural_accounting.py" --root "$D/structural_fallback_scaling_v1" --stage complete --output "$D/structural_fallback_accounting_v1.json"
"$A/.venv/bin/python" "$E/assess_fallback_scaling.py" --root "$D/structural_fallback_scaling_v1" --accounting "$D/structural_fallback_accounting_v1.json" --output "$D/structural_fallback_assessment_v1"
