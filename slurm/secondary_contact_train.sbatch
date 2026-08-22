#!/usr/bin/env bash
#SBATCH --job-name=sc-train
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
if [[ -n "${CHUNKS_OVERRIDE:-}" ]]; then
  CHUNKS="$CHUNKS_OVERRIDE"
elif [[ "${DATASET:-chromosome}" == "mixed_sfs" ]]; then
  CHUNKS='data/processed/secondary_contact/sfs_chunks/chunk_*.npz,data/processed/secondary_contact/focused_sfs_chunks/chunk_*.npz'
elif [[ "${DATASET:-chromosome}" == "focused_sfs" ]]; then
  CHUNKS='data/processed/secondary_contact/focused_sfs_chunks/chunk_*.npz'
elif [[ "${DATASET:-chromosome}" == "arg" ]]; then
  CHUNKS='data/processed/secondary_contact/arg_training.npz'
elif [[ "${DATASET:-chromosome}" == "sfs" ]]; then
  CHUNKS='data/processed/secondary_contact/sfs_chunks/chunk_*.npz'
else
  CHUNKS='data/processed/secondary_contact/chunks/chunk_*.npz'
fi
EXTRA_ARGS=()
if [[ -n "${DOMAIN_FILTER:-}" ]]; then
  EXTRA_ARGS+=(--domain-filter "$DOMAIN_FILTER")
fi
if [[ -n "${NUISANCE_FILTER:-}" ]]; then
  EXTRA_ARGS+=(--nuisance-filter "$NUISANCE_FILTER")
fi
.venv/bin/python scripts/train_secondary_contact.py \
  --chunks "$CHUNKS" \
  --out "artifacts/task2_secondary_contact/model_${DATASET:-chromosome}_${BACKEND:-transformer}_${SEED:-7701}.pt" \
  --context-backend "${BACKEND:-transformer}" \
  --seed "${SEED:-7701}" \
  --epochs "${EPOCHS:-120}" \
  --patience "${PATIENCE:-20}" \
  --width "${WIDTH:-128}" \
  --score-weight "${SCORE_WEIGHT:-0.15}" \
  --objective "${OBJECTIVE:-posterior}" \
  --spectrum-transform "${SPECTRUM_TRANSFORM:-log_l2}" \
  "${EXTRA_ARGS[@]}"
