#!/usr/bin/env bash
#SBATCH --job-name=task11-development
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-19%4
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v4"
export PYTHONPATH="$PWD:$PWD/src:$PWD/scripts"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if (( SLURM_ARRAY_TASK_ID < 16 )); then split=train; index=$SLURM_ARRAY_TASK_ID; else split=calibration; index=$((SLURM_ARRAY_TASK_ID - 16)); fi
cohort=$(printf '%s_%02d' "$split" "$index")
"$PROJECT_ROOT/.venv/bin/python" scripts/develop_task11_pedigree.py run \
 --protocol "$TASK_ROOT/protocol_v4.json" --split "$split" --index "$index" \
 --map "$PROJECT_ROOT/external/dog-recombination/campbell2016_average.plink.map" \
 --hapibd "$PROJECT_ROOT/external/hap-ibd/hap-ibd.jar" \
 --out "$TASK_ROOT/results_v4/$cohort" --baseline "$TASK_ROOT/results_v3/$cohort"
