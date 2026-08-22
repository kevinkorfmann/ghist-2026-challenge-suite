#!/usr/bin/env bash
#SBATCH --job-name=bn-observed
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
VCF="data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz"

.venv/bin/python -m ghist_bottleneck.cli extract-blocks \
  --vcf "$VCF" \
  --out data/processed/task1_bottleneck/observed.blocks.npz \
  --sequence-length 100000000 \
  --block-length 2000000

.venv/bin/python -m ghist_bottleneck.cli sfs-reference \
  --vcf "$VCF" \
  --out artifacts/task1_bottleneck/moments_unfolded.json \
  --restarts 64

.venv/bin/python -m ghist_bottleneck.cli sfs-reference \
  --vcf "$VCF" \
  --out artifacts/task1_bottleneck/moments_folded.json \
  --restarts 64 \
  --folded
