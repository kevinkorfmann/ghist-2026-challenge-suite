#!/usr/bin/env bash
#SBATCH --job-name=bn-mld-prep
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
.venv/bin/python -m ghist_bottleneck.momentsld_reference prepare \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz
