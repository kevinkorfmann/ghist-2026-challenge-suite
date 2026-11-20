#!/bin/bash
#SBATCH --job-name=task10-multiway
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/multiway-train-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/multiway-source"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$S/experiments/task10_development_2026_11_20/multiway_candidate.py" --asset-root "$A" --source-root "$S" --phase-model "$D/phase_training/phase_model.joblib" --protocol "$D/multiway_protocol.json" --out-dir "$D/multiway_training"
