#!/bin/bash
#SBATCH --job-name=task9-current-labels
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/current-%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/current-%j.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/scientific-relatedness-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
"$TASK_ROOT/.venv/bin/python" scripts/train_task9_current_candidate.py \
  --baseline-model "$TASK_ROOT/artifacts/task9_ibd_segments/model_v3_shape/task9_model.joblib" \
  --calibration "$TASK_ROOT"/artifacts/task9_ibd_segments/calibration/replicate_{0,1,2,3,4,5}.seg \
  --sex-map "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" \
  --out-dir artifacts/task9_current_contract --workers 4
for index in {0..5}; do
  if (( index < 3 )); then scenario=interference; else scenario=poisson; fi
  seed=$((2026112090 + index))
  "$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py \
    --root "$TASK_ROOT" --segments "artifacts/task9_${scenario}_${index}.seg" \
    --kind task9 --scenario "$scenario" --seed "$seed" \
    --task9-model artifacts/task9_current_contract/model.joblib \
    --out "results/task9_current_${scenario}_${index}.json"
  if (( index >= 3 )); then
    "$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py \
      --root "$TASK_ROOT" --segments "artifacts/task9_${scenario}_${index}.seg" \
      --kind task9 --scenario "$scenario" --seed "$seed" \
      --out "results/task9_${scenario}_${index}.json"
  fi
done
