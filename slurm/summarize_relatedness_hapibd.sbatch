#!/usr/bin/env bash
#SBATCH --job-name=relatedness-ibd-summary
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

case "$SLURM_ARRAY_TASK_ID" in
  0)
    ROOT="artifacts/task10_large_family/hapibd"
    MAP="external/ped-sim-map/refined_average.plink.map"
    ;;
  1)
    ROOT="artifacts/task11_multifamily_inbreeding/hapibd"
    MAP="external/dog-recombination/campbell2016_average.plink.map"
    ;;
  2)
    ROOT="artifacts/task12_large_family_background/hapibd"
    MAP="external/ped-sim-map/refined_average.plink.map"
    ;;
esac

PYTHONPATH=src .venv/bin/python scripts/summarize_hapibd.py \
  --ibd "$ROOT/segments.ibd.gz" \
  --hbd "$ROOT/segments.hbd.gz" \
  --map "$MAP" \
  --out-dir "$ROOT/summary"
