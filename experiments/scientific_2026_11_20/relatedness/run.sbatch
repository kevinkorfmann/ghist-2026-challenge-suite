#!/bin/bash
#SBATCH --job-name=ghist-relatedness-science
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-11%4
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
SCIENCE_ROOT="$TASK_ROOT/reviews/scientific-relatedness-20261120"
cd "$SCIENCE_ROOT"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
index=$SLURM_ARRAY_TASK_ID
if (( index < 6 )); then
  kind=task9
  definition=configs/task9_pedigrees.def
  map="$TASK_ROOT/external/ped-sim-map/refined_mf.simmap"
  if (( index < 3 )); then
    scenario=interference
    flags=(--intf "$TASK_ROOT/external/ped-sim/interfere/nu_p_campbell.tsv")
  else
    scenario=poisson
    flags=(--pois)
  fi
elif (( index < 9 )); then
  kind=human
  scenario=poisson
  definition=configs/task10_human_pedigrees.def
  map="$TASK_ROOT/external/ped-sim-map/refined_mf.simmap"
  flags=(--pois)
else
  kind=dog
  scenario=poisson
  definition=configs/task11_pedigrees.def
  map="$TASK_ROOT/external/dog-recombination/campbell2016_average.simmap"
  flags=(--pois)
fi
seed=$((2026112090 + index))
"$TASK_ROOT/external/ped-sim/ped-sim" -d "$definition" -m "$map" -o "artifacts/${kind}_${scenario}_${index}" "${flags[@]}" --seed "$seed"
"$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py --root "$TASK_ROOT" --segments "artifacts/${kind}_${scenario}_${index}.seg" --kind "$kind" --scenario "$scenario" --seed "$seed" --out "results/${kind}_${scenario}_${index}.json"
