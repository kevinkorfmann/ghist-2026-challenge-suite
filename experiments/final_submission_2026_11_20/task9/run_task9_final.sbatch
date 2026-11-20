#!/bin/bash
#SBATCH --job-name=ghist-task9-final
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task9-root-20261120-1914/job-%j.log
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
TASK_PYTHON="$TASK_ROOT/.venv/bin/python"
TASK_PACKAGE="$TASK_ROOT/reviews/task9-offline-package-v2-20261120/task9-package-v2"
TASK_FINAL="$TASK_ROOT/final_runs/task9-root-20261120-1914"
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_FINAL"
"$TASK_PYTHON" "$TASK_PACKAGE/scripts/prepare_task9_submission_package.py" run \
  --package "$TASK_PACKAGE" --artifact-root "$TASK_ROOT" \
  --sex-map "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" \
  --observed "$TASK_FINAL/relatedness1.final.csv" \
  --input-provenance "$TASK_FINAL/input_provenance.json" \
  --run-dir "$TASK_FINAL/run_v1" --workers 2 \
  --expected-package-sha256 4c0ec8f9ae12b00b77506b803812c3378f1d43bf873ff206a1db2bffeca5da9e
