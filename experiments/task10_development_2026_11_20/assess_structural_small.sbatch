#!/bin/bash
#SBATCH --job-name=task10-small-assess
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/small-assess-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/structural-qualification-source-v3"
export PYTHONPATH="$Q:$S/experiments/task10_development_2026_11_20:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$D/collect_structural_accounting.py" --root "$D/structural_qualification_v3" --stage small --output "$D/structural_small_accounting_v3.json"
"$A/.venv/bin/python" "$Q/assess_structural_qualification.py" --root "$D/structural_qualification_v3" --stage small --protocol "$Q/structural_qualification_protocol_v3.json" --accounting "$D/structural_small_accounting_v3.json" --output "$D/structural_small_assessment_v3"
