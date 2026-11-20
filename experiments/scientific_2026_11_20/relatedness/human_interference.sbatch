#!/bin/bash
#SBATCH --job-name=human-independent-interference
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=12-14%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/human-%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/human-%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/scientific-relatedness-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
index=$SLURM_ARRAY_TASK_ID
seed=$((2026112090 + index))
"$TASK_ROOT/external/ped-sim/ped-sim" -d configs/task10_human_pedigrees.def -m "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" -o "artifacts/human_interference_${index}" --intf "$TASK_ROOT/external/ped-sim/interfere/nu_p_campbell.tsv" --seed "$seed"
"$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py --root "$TASK_ROOT" --segments "artifacts/human_interference_${index}.seg" --kind human --scenario interference --seed "$seed" --out "results/human_interference_${index}.json"
