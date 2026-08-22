#!/usr/bin/env bash
#SBATCH --job-name=growth-bgs-resources
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

PYTHONPATH=src .venv/bin/python scripts/audit_growth_bgs.py \
  --vcf data/raw/demography/growth_background_selection/testing/growth.testing.vcf.gz \
  --coding-bed data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed \
  --recombination-map data/raw/shared/demography/comeron2012v2_combined_map_dm6.txt \
  --out artifacts/task3_growth_bgs/observed_resource_audit.json
