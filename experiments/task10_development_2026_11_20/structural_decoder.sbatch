#!/bin/bash
#SBATCH --job-name=task10-nj10-diagnostic
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/nj10-diagnostic-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/structural-decoder-source-v1"
P="$S/experiments/task10_development_2026_11_20"
export PYTHONPATH="$Q:$P:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$Q/diagnose_structural_decoder.py" --source-root "$S" --protocol "$Q/structural_decoder_protocol.json" --panel-root "$D/detected_panel_v2" --output "$D/structural_decoder_v1"
