#!/usr/bin/env bash
#SBATCH --job-name=splitmig-block-extract
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/extract_split_migration_bgs_blocks.py \
  --vcf data/raw/demography/split_migration_background_selection/testing/split_migration.testing.vcf.gz \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --block-size 1000000 \
  --out artifacts/task4_split_migration_bgs/blocks_1mb.npz \
  --summary artifacts/task4_split_migration_bgs/blocks_1mb_summary.json

