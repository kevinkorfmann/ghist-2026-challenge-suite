#!/usr/bin/env bash
#SBATCH --job-name=sc-gen-bootstrap-aggregate
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/summarize_secondary_contact_generator_bootstrap.py \
  --inputs "${INPUTS:-artifacts/task2_secondary_contact/generator_bootstrap/replicate_*.json}" \
  --out "${OUT:-artifacts/task2_secondary_contact/generator_bootstrap_summary.json}" \
  --observed \
  "${MAINLAND:-20127.674}" "${ISLAND:-79105.691}" \
  "${SPLIT:-4346.354}" "${CONTACT:-1395.763}" "${MIGRATION:-0.000207240}"
