#!/bin/bash
#SBATCH --job-name=ghist-task14-pragmatic-final
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:55:00
# Root supplies --output to sbatch and the five required environment values below.
source /etc/profile
set -euo pipefail
: "${ROOT_DECISION:?root-owned decision JSON is required}"
: "${DECISION_SHA:?verified root decision SHA256 is required}"
: "${FINAL_INPUT_DIRECTORY:?root-verified BAM/popfile directory is required}"
: "${FINAL_REFERENCE:?root-verified matching current reference is required}"
: "${FINAL_OUTPUT:?new root-owned output directory is required}"
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
export PYTHONPATH="$TASK_ROOT/source_v12/src:$TASK_ROOT/source_v12"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v12"
"$TASK_PYTHON" scripts/lowcoverage_predict_pragmatic.py \
 --frozen-root "$TASK_ROOT" \
 --freeze-sha 7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb \
 --root-decision "$ROOT_DECISION" --decision-sha "$DECISION_SHA" \
 --directory "$FINAL_INPUT_DIRECTORY" --reference "$FINAL_REFERENCE" --output "$FINAL_OUTPUT"
