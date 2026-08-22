#!/usr/bin/env bash
#SBATCH --job-name=growth-bgs-summary
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/summarize_growth_bgs_audit.py \
  --generic-audit artifacts/task3_growth_bgs/observed_vcf_audit.json \
  --resource-audit artifacts/task3_growth_bgs/observed_resource_audit.json \
  --out artifacts/task3_growth_bgs/observed_audit_summary.json
