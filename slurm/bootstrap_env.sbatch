#!/usr/bin/env bash
#SBATCH --job-name=ghist-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
python -m venv .venv
.venv/bin/pip install --upgrade pip wheel setuptools
.venv/bin/pip install --extra-index-url https://download.pytorch.org/whl/cu128 -e .
.venv/bin/python -c 'import torch, msprime, pysam; print(torch.__version__, torch.cuda.is_available())'

