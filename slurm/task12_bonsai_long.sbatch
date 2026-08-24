#!/usr/bin/env bash
#SBATCH --job-name=task12-bonsai-long
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/reconstruct_bonsai_pedigree.py \
  --ibd artifacts/task12_large_family_background/hapibd/segments.ibd.gz \
  --map external/ped-sim-map/refined_average.plink.map \
  --threshold-cm 40 \
  --out artifacts/task12_large_family_background/hapibd/bonsai_long/summary.json
