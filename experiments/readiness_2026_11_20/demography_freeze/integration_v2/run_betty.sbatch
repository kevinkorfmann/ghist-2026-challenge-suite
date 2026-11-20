#!/bin/bash
#SBATCH --job-name=ghist-freeze-checks
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
GHIST_REVIEW_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
GHIST_REVIEW_CODE="$GHIST_REVIEW_ROOT/reviews/demography-integrity-v2-integration-20261120"
export PYTHONPATH="$GHIST_REVIEW_ROOT/reviews/demography-integrity-integration-20261120/.lint-deps:$GHIST_REVIEW_ROOT/reviews/2026-11-20/.test-deps"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$GHIST_REVIEW_CODE"
"$GHIST_REVIEW_ROOT/.venv/bin/python" run_checks.py "$GHIST_REVIEW_CODE" betty_checks
