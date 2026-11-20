#!/usr/bin/env bash
#SBATCH --job-name=soft-sweep-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/audit_single_sweep.py \
  --vcf data/raw/sweeps/soft_sweep/testing/softsweep.testing.vcf.gz \
  --accessible-bed data/raw/shared/sweeps/testing.accessible.bed \
  --score-field soft_composite_score \
  --out-dir artifacts/task6_soft_sweep/baseline_scans \
  --summary artifacts/task6_soft_sweep/audit_summary.json
