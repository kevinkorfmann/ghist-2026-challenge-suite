#!/usr/bin/env bash
#SBATCH --job-name=human-pair-topology-benchmark
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/benchmark_bonsai_pair_topology.py \
  --segments artifacts/task10_human_calibration/replicate_6.seg \
  --threshold-cm 5 \
  --pairs-per-label 100 \
  --out artifacts/task10_human_calibration/bonsai_pair_benchmark/summary.json
