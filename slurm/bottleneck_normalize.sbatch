#!/usr/bin/env bash
#SBATCH --job-name=bn-vcf-normalize
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python -m ghist_bottleneck.vcf_normalize \
  --source data/raw/task1_bottleneck/bottleneck.testing.vcf.gz \
  --destination data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --report artifacts/task1_bottleneck/vcf_normalization.json
