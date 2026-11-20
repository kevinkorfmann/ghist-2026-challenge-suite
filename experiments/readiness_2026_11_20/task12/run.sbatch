#!/bin/bash
#SBATCH --job-name=task12-reference-refit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task12-20261120/job-%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task12-20261120/job-%j.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/readiness-task12-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1
PROTOCOL=experiments/readiness_2026_11_20/task12/protocol.json
PYTHON="$TASK_ROOT/.venv/bin/python"
# The candidate is fitted before any fresh holdout is generated.
"$PYTHON" scripts/train_task12_reference_candidate.py --root "$TASK_ROOT" \
  --protocol "$PROTOCOL" --out-dir artifacts/task12_reference_v2
for SEED in 2026114011 2026114012 2026114013 2026114014; do
  test ! -e "artifacts/holdout_${SEED}.seg"
  "$TASK_ROOT/external/ped-sim/ped-sim" -d configs/task10_human_pedigrees.def \
    -m "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" \
    --intf "$TASK_ROOT/external/ped-sim/interfere/nu_p_campbell.tsv" \
    --seed "$SEED" -o "artifacts/holdout_${SEED}"
done
"$PYTHON" scripts/validate_task12_reference_candidate.py --root "$TASK_ROOT" \
  --protocol "$PROTOCOL" --candidate-dir artifacts/task12_reference_v2 \
  --holdouts artifacts/holdout_2026114011.seg artifacts/holdout_2026114012.seg \
             artifacts/holdout_2026114013.seg artifacts/holdout_2026114014.seg \
  --out-dir results
