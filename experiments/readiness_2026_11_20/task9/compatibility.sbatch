#!/usr/bin/env bash
#SBATCH --job-name=task9-runtime-compatibility
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task9-compat-20261120/job-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
SCIENTIFIC_ROOT="$TASK_ROOT/reviews/scientific-relatedness-20261120"
cd "$TASK_ROOT/reviews/readiness-task9-compat-20261120"
export PYTHONPATH="$PWD/src:$PWD" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
"$TASK_ROOT/.venv/bin/python" scripts/validate_task9_operationally.py \
  --release configs/task9_current_release.json --artifact-root "$TASK_ROOT" \
  --model "$SCIENTIFIC_ROOT/artifacts/task9_current_contract/model.joblib" \
  --sex-map "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" \
  --segments "$SCIENTIFIC_ROOT/artifacts/task9_confirmation_interference.seg" \
             "$SCIENTIFIC_ROOT/artifacts/task9_confirmation_poisson.seg" \
  --out-dir results/versioned
