#!/bin/bash
#SBATCH --job-name=task10-detected-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/detected-fit-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
P="$S/experiments/task10_development_2026_11_20"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$P/detected_family_panel.py" fit --asset-root "$S/assets" --source-root "$S" --phase-model "$S/models/phase_model.joblib" --protocol "$P/detected_panel_protocol_v2.json" --out-root "$D/detected_panel_v2"
