#!/usr/bin/env bash
#SBATCH --job-name=sc-generator-check
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/validate_secondary_contact_generator_spectrum.py \
  --out artifacts/task2_secondary_contact/generator_spectrum_validation.json \
  --replicates "${REPLICATES:-20}" \
  --sequence-length "${SEQUENCE_LENGTH:-500000}" \
  --seed "${SEED:-20260907}"
