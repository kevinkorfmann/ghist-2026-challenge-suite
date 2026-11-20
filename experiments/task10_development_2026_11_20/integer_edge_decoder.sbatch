#!/bin/bash
#SBATCH --job-name=task10-integer-edges
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/integer-edges-%j.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/detected-panel-source"
Q="$D/integer-edge-source-v1"
P="$S/experiments/task10_development_2026_11_20"
export PYTHONPATH="$Q:$P:$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$A/.venv/bin/python" "$Q/integer_edge_decoder.py" --source-root "$S" --protocol "$Q/integer_edge_protocol.json" --panel-root "$D/detected_panel_v2" --initial-tree-root "$D/forest_distance_decoder_v1" --output "$D/integer_edge_decoder_v1"
