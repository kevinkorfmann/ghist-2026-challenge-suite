#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python scripts/audit_single_sweep.py \
  --vcf data/raw/sweeps/single_sweep/testing/onesweep.testing.vcf.gz \
  --accessible-bed data/raw/shared/sweeps/testing.accessible.bed \
  --out-dir artifacts/task5_single_sweep/baseline_scans \
  --summary artifacts/task5_single_sweep/audit_summary.json

