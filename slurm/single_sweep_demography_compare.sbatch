#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-demography
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/compare_single_sweep_demographies.py \
  --vcf data/raw/sweeps/single_sweep/testing/onesweep.testing.vcf.gz \
  --accessible-bed data/raw/shared/sweeps/testing.accessible.bed \
  --out artifacts/task5_single_sweep/demography_comparison.json

