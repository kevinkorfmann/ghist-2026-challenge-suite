#!/bin/bash
#SBATCH --job-name=task9-independent-rescore
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --array=0-5%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/rescore-%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/scientific-relatedness-20261120/logs/rescore-%A_%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/reviews/scientific-relatedness-20261120"
export PYTHONPATH="$PWD/src:$PWD:$PWD/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
index=$SLURM_ARRAY_TASK_ID
if (( index < 3 )); then scenario=interference; else scenario=poisson; fi
"$TASK_ROOT/.venv/bin/python" scripts/validate_relatedness_scientifically.py --root "$TASK_ROOT" --segments "artifacts/task9_${scenario}_${index}.seg" --kind task9 --scenario "$scenario" --seed "$((2026112090 + index))" --out "results/task9_${scenario}_${index}.json"
