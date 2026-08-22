#!/usr/bin/env bash
#SBATCH --job-name=sc-tabular
#SBATCH --partition=genoa-lrg-mem
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=15872M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
.venv/bin/python scripts/train_secondary_contact_tabular.py \
  --chunks "${CHUNKS:-data/processed/secondary_contact/generator_sfs_chunks/chunk_*.npz}" \
  --out "${OUT:-artifacts/task2_secondary_contact/model_generator_sfs_tabular_0.joblib}" \
  --domain-filter "${DOMAIN_FILTER:-0}" \
  --pca-components "${PCA_COMPONENTS:-96}" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --seed "${SEED:-20260923}"
