#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/%A_%a.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v5"
export PYTHONPATH="$PWD:$PWD/src:$PWD/scripts"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cohort=$(printf 'confirmation_%02d' "$SLURM_ARRAY_TASK_ID")
"$PROJECT_ROOT/.venv/bin/python" scripts/develop_task11_pedigree.py run --protocol "$TASK_ROOT/protocol_v5.json" --split confirmation --index "$SLURM_ARRAY_TASK_ID" --candidate "$TASK_ROOT/candidate_v5/candidate.json" --map "$PROJECT_ROOT/external/dog-recombination/campbell2016_average.plink.map" --hapibd "$PROJECT_ROOT/external/hap-ibd/hap-ibd.jar" --out "$TASK_ROOT/results_v5/$cohort"
