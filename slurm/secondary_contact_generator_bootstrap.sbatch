#!/usr/bin/env bash
#SBATCH --job-name=sc-gen-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
OUTPUT_ROOT="${OUTPUT_ROOT:-artifacts/task2_secondary_contact/generator_bootstrap}"
.venv/bin/python scripts/bootstrap_secondary_contact_generator_fit.py \
  --out "$OUTPUT_ROOT/replicate_${SLURM_ARRAY_TASK_ID}.json" \
  --replicate "$SLURM_ARRAY_TASK_ID" \
  --truth \
  "${ANCESTRAL:-32865}" "${MAINLAND:-20100}" "${ISLAND:-79400}" \
  "${SPLIT:-4350}" "${CONTACT:-1390}" "${MIGRATION:-0.000208}" \
  --blocks "${BLOCKS:-10}" \
  --block-length "${BLOCK_LENGTH:-10000000}" \
  --projection ${PROJECTION:-40 30} \
  --restarts "${RESTARTS:-12}" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --maxiter "${MAXITER:-400}" \
  --seed "${SEED:-20260919}"
