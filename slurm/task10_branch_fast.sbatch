#!/usr/bin/env bash
#SBATCH --job-name=task10-branch-fast
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
root="artifacts/task10_large_family/hapibd"

PYTHONPATH=src:scripts .venv/bin/python scripts/reconstruct_human_branch_pedigrees.py \
  --ibd "$root/segments.ibd.gz" \
  --map external/ped-sim-map/refined_average.plink.map \
  --clusters "$root/branch_clusters/summary.json" \
  --threshold-cm 10 \
  --max-add-degree 3 \
  --max-rel-types 1 \
  --out "$root/branch_bonsai_fast/summary.json"
