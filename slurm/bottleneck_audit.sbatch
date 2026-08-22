#!/usr/bin/env bash
#SBATCH --job-name=bn-vcf-audit
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

.venv/bin/python scripts/audit_vcf.py \
  --vcf data/raw/task1_bottleneck/bottleneck.testing.vcf.gz \
  --out artifacts/task1_bottleneck/observed_vcf_audit.json
