#!/usr/bin/env bash
#SBATCH --job-name=task11-complete-smoke
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
#SBATCH --array=0-1
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-complete-20261120/logs/%A_%a.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-complete-20261120/logs/%A_%a.err
set -euo pipefail
module load python/3.11.11
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
TASK_ROOT="$PROJECT_ROOT/experiments/task11-complete-20261120"
cd "$TASK_ROOT/source"
export PYTHONPATH="$PWD:$PWD/src:$PWD/scripts"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
seeds=(2026112711 2026112811)
"$PROJECT_ROOT/.venv/bin/python" scripts/validate_task11_complete_pipeline.py \
  --map "$PROJECT_ROOT/external/dog-recombination/campbell2016_average.plink.map" \
  --hapibd "$PROJECT_ROOT/external/hap-ibd/hap-ibd.jar" \
  --model "$PROJECT_ROOT/artifacts/task11_multifamily_inbreeding/pair_model_hbd5/model.joblib" \
  --seed "${seeds[$SLURM_ARRAY_TASK_ID]}" \
  --out "$TASK_ROOT/results/seed_${seeds[$SLURM_ARRAY_TASK_ID]}"
