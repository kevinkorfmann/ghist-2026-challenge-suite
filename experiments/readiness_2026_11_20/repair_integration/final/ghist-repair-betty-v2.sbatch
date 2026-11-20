#!/usr/bin/env bash
#SBATCH --job-name=ghist-repair-tests
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/repair-integration-20261120/slurm-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/repair-integration-20261120"
export PYTHONPATH="$PWD/src:$PWD:$TASK_ROOT/reviews/2026-11-20/.test-deps"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/.mplconfig"
"$TASK_ROOT/.venv/bin/python" -m pytest -q > betty_pytest_v2.log 2>&1

PYTHONPATH="$TASK_ROOT/reviews/2026-11-20/.test-deps" "$TASK_ROOT/.venv/bin/python" ghist_repair_evidence.py "$PWD" "$PWD/source_snapshot.json" > betty_environment_v2.json
