#!/bin/bash
#SBATCH --job-name=lc13-prag-final
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:55:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/task13-prag-final-%j.log
source /etc/profile >/dev/null
set -euo pipefail
: "${ROOT_DECISION:?}" "${DECISION_SHA:?}" "${FINAL_INPUT_DIRECTORY:?}" "${FINAL_REFERENCE:?}" "${FINAL_OUTPUT:?}"
cd /vast/projects/smathi/cohort/korfmann/ghist2026
export PYTHONPATH="$PWD/experiments/lowcoverage-20261120/source_v12/src:$PWD/experiments/lowcoverage-20261120/source_v12"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python experiments/lowcoverage-20261120/lowcoverage_predict_task13_pragmatic.py \
  --frozen-root "$PWD/experiments/lowcoverage-20261120" \
  --root-decision "$ROOT_DECISION" --decision-sha "$DECISION_SHA" \
  --directory "$FINAL_INPUT_DIRECTORY" --reference "$FINAL_REFERENCE" --output "$FINAL_OUTPUT"
