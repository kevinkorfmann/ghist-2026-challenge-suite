#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/pragmatic_corrected_inference_%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/pragmatic_corrected_inference_%j.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
ROLE="$3"
case "$ROLE" in base|context) ;; *) echo 'Unknown fixed model role' >&2; exit 2 ;; esac
cd "$TASK_ROOT/source_v5"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[[ "$(sha256sum "$TASK_ROOT/pragmatic_corrected_bindings_v1.py" | cut -d' ' -f1)" == 580c492cd8997821669ee4677df11a0123c21e79ad95b292704ecaae011d3c9f ]]
[[ "$(sha256sum "$TASK_ROOT/prepare_pragmatic_corrected_v1.py" | cut -d' ' -f1)" == 69b288845b8a44f2fc406f43b90daa6b3a05ca3f9570c7b0645b9c238b95f632 ]]
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/prepare_pragmatic_corrected_v1.py" "$TASK_ROOT" \
  --vcf "$1" --out "$2" --candidate "$TASK_ROOT/pragmatic_corrected_${ROLE}_candidate_v1.json" \
  --candidate-sha256 "$4" --runner-sha256 580c492cd8997821669ee4677df11a0123c21e79ad95b292704ecaae011d3c9f --dataset final
