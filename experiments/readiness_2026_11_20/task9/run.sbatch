#!/bin/bash
#SBATCH --job-name=task9-production-readiness
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:35:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task9-20261120/job-%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task9-20261120/job-%j.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
SCIENTIFIC_ROOT="$TASK_ROOT/reviews/scientific-relatedness-20261120"
cd "$TASK_ROOT/reviews/readiness-task9-20261120"
export PYTHONPATH="$PWD/src:$PWD"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
MODEL="$SCIENTIFIC_ROOT/artifacts/task9_current_contract/model.joblib"
MAP="$TASK_ROOT/external/ped-sim-map/refined_mf.simmap"
SEGMENTS=("$SCIENTIFIC_ROOT/artifacts/task9_confirmation_interference.seg" "$SCIENTIFIC_ROOT/artifacts/task9_confirmation_poisson.seg")
"$TASK_ROOT/.venv/bin/python" scripts/validate_task9_operationally.py \
  --model "$MODEL" --sex-map "$MAP" --artifact-root "$TASK_ROOT" \
  --segments "${SEGMENTS[@]}" --out-dir results/unversioned
"$TASK_ROOT/.venv/bin/python" scripts/prepare_task9_release.py \
  --model "$MODEL" --training "$SCIENTIFIC_ROOT/artifacts/task9_current_contract/training.json" \
  --confirmation "$SCIENTIFIC_ROOT/results/task9_confirmation_current_interference.json" \
                 "$SCIENTIFIC_ROOT/results/task9_confirmation_current_poisson.json" \
  --operational-report results/unversioned/operational.json --sex-map "$MAP" \
  --artifact-root "$TASK_ROOT" --manifest configs/task9_current_release.json
"$TASK_ROOT/.venv/bin/python" scripts/validate_task9_operationally.py \
  --release configs/task9_current_release.json --artifact-root "$TASK_ROOT" \
  --model "$MODEL" --sex-map "$MAP" --segments "${SEGMENTS[@]}" --out-dir results/versioned
