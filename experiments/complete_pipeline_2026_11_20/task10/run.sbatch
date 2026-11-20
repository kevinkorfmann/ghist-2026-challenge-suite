#!/usr/bin/env bash
#SBATCH --job-name=task10-complete-20261120
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=03:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/artifacts/task10-complete-20261120/slurm-%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/artifacts/task10-complete-20261120/slurm-%j.err
set -euo pipefail
module load python/3.11.11
ASSETS=/vast/projects/smathi/cohort/korfmann/ghist2026
SNAPSHOT="$ASSETS/artifacts/task10-complete-20261120/source"
export PYTHONPATH="$SNAPSHOT/src:$SNAPSHOT/scripts"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd "$SNAPSHOT"
"$ASSETS/.venv/bin/python" scripts/validate_task10_complete_pipeline.py \
  --project-root "$SNAPSHOT" --asset-root "$ASSETS" \
  --out-dir "$ASSETS/artifacts/task10-complete-20261120/seed_2026112010" \
  --seed 2026112010 --workers 4
