#!/usr/bin/env bash
#SBATCH --job-name=sc-arg-fit
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
test -n "${FEATURES:-}"
test -n "${OUT:-}"
.venv/bin/python scripts/fit_secondary_contact_arg.py \
  --features "$FEATURES" \
  --out "$OUT" \
  --initial \
  "${ANCESTRAL:-16578}" "${MAINLAND:-18725}" "${ISLAND:-51450}" \
  "${SPLIT:-4638}" "${CONTACT:-1183}" "${MIGRATION:-0.000245}" \
  --seed "${SEED:-20260822}" \
  ${LOCAL_ONLY:+--local-only}
