#!/usr/bin/env bash
#SBATCH --job-name=bn-mom-chrom
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
ENGINE="${ENGINE:-moments}"
OBSERVED_FIT="${OBSERVED_FIT:-artifacts/task1_bottleneck/moments_folded.json}"
OUT="${OUT:-artifacts/task1_bottleneck/moments_chromosome_validation.json}"
.venv/bin/python scripts/evaluate_moments_chromosomes.py \
  --dataset data/simulations/bottleneck/chromosome_validation.npz \
  --split-record artifacts/task1_bottleneck/chromosome_validation_7732378.json \
  --observed-fit "$OBSERVED_FIT" \
  --out "$OUT" \
  --workers 16 \
  --restarts 12 \
  --engine "$ENGINE"
