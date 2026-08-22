#!/usr/bin/env bash
#SBATCH --job-name=bn-genotype-matrix
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src .venv/bin/python -m ghist_bottleneck.genotype_matrix \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --output-dir data/processed/task1_bottleneck/fastcxt_input \
  --expected-records 1194350
