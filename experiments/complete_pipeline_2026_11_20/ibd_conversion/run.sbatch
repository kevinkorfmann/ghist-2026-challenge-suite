#!/usr/bin/env bash
#SBATCH --job-name=ibd-conversion-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/ibd-conversion-20261120/logs/%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/ibd-conversion-20261120/logs/%j.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
TASK_ROOT="$PROJECT_ROOT/experiments/ibd-conversion-20261120"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT_ROOT/.venv/bin/python" "$TASK_ROOT/audit.py" \
  --historical-root "$PROJECT_ROOT/experiments/task11-complete-20261120" \
  --corrected "$TASK_ROOT/corrected_converter.py" \
  --out "$TASK_ROOT/results"
