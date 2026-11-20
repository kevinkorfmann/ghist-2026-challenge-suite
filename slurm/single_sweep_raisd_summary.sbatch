#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-raisd-summary
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
RESULT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/raisd"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/summarize_single_sweep_raisd.py \
  --report "$RESULT_ROOT/RAiSD_Report.task5_w50" \
  --info "$RESULT_ROOT/RAiSD_Info.task5_w50" \
  --site-report "$RESULT_ROOT/RAiSD_SiteReport.task5_w50" \
  --job-id 7770050 \
  --out artifacts/task5_single_sweep/raisd_summary.json

