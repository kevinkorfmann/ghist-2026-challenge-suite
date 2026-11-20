#!/bin/bash
#SBATCH --job-name=task10-inference
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:30:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/final-inference-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/structural-qualification-source-v3"
H="$D/structural-inference-handoff-v2"
export PYTHONPATH="$Q:$S/experiments/task10_development_2026_11_20:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
# Root supplies input path, SHA256, isolated output and optional explicit missing policy.
"$A/.venv/bin/python" "$H/structural_inference_handoff_v2.py" --handoff-source "$H" "$@"
