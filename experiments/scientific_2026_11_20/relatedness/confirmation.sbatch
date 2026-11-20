#!/bin/bash
#SBATCH --job-name=task9-untouched-confirmation
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=0-1
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/confirmation-%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/confirmation-%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/scientific-relatedness-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
index=$SLURM_ARRAY_TASK_ID
seed=$((2026112307 + index))
if (( index == 0 )); then
  scenario=interference
  flags=(--intf "$TASK_ROOT/external/ped-sim/interfere/nu_p_campbell.tsv")
else
  scenario=poisson
  flags=(--pois)
fi
prefix="artifacts/task9_confirmation_${scenario}"
"$TASK_ROOT/external/ped-sim/ped-sim" -d configs/task9_pedigrees.def -m "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" -o "$prefix" "${flags[@]}" --seed "$seed"
for version in frozen current; do
  flags=()
  if [[ "$version" == current ]]; then flags=(--task9-model artifacts/task9_current_contract/model.joblib); fi
  "$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py \
    --root "$TASK_ROOT" --segments "${prefix}.seg" --kind task9 \
    --scenario "confirmation_${scenario}" --seed "$seed" "${flags[@]}" \
    --out "results/task9_confirmation_${version}_${scenario}.json"
done
