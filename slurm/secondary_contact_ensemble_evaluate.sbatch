#!/usr/bin/env bash
#SBATCH --job-name=sc-ensemble-evaluate
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
.venv/bin/python scripts/evaluate_secondary_contact_ensemble.py \
  --evaluations "${EVALUATIONS:-artifacts/task2_secondary_contact/model_generator_sfs_transformer_*.evaluation.npz}" \
  --out "${OUT:-artifacts/task2_secondary_contact/generator_sfs_ensemble_evaluation.json}" \
  --domain-filter "${DOMAIN_FILTER:-0}" \
  --relative-selection-cutoff "${RELATIVE_SELECTION_CUTOFF:-1.25}" \
  --score-power "${SCORE_POWER:-2}"
