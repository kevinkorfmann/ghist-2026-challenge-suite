#!/usr/bin/env bash
#SBATCH --job-name=human-branch-pedigrees
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=05:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  root="artifacts/task10_large_family/hapibd"
  threshold_cm=10
else
  root="artifacts/task12_large_family_background/hapibd"
  threshold_cm=20
fi

PYTHONPATH=src:scripts .venv/bin/python scripts/cluster_human_family_branches.py \
  --pairs "$root/summary/pair_features.csv" \
  --threshold 0.12 \
  --out "$root/branch_clusters/summary.json"

PYTHONPATH=src:scripts .venv/bin/python scripts/reconstruct_human_branch_pedigrees.py \
  --ibd "$root/segments.ibd.gz" \
  --map external/ped-sim-map/refined_average.plink.map \
  --clusters "$root/branch_clusters/summary.json" \
  --threshold-cm "$threshold_cm" \
  --out "$root/branch_bonsai/summary.json"
