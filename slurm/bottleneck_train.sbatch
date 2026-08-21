#!/usr/bin/env bash
#SBATCH --job-name=bn-fit
#SBATCH --partition=dgx-b200
#SBATCH --gres=gpu:B200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
module load python/3.11.11 cuda/12.8.1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/merge_simulations.py \
  --input 'data/simulations/bottleneck/part_*.npz' \
  --output data/simulations/bottleneck/all.npz
.venv/bin/python -m ghist_bottleneck.cli train \
  --dataset data/simulations/bottleneck/all.npz \
  --out-dir checkpoints/bottleneck

