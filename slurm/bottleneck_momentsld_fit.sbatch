#!/usr/bin/env bash
#SBATCH --job-name=bn-mld-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
MODE="${MODE:-joint}"
EXTRA=()
if [[ "$MODE" == "ld_only" ]]; then
  EXTRA+=(--ld-only)
elif [[ "$MODE" != "joint" ]]; then
  echo "MODE must be joint or ld_only" >&2
  exit 2
fi
.venv/bin/python -m ghist_bottleneck.momentsld_reference fit \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --region-dir artifacts/task1_bottleneck/momentsld_regions \
  --out "artifacts/task1_bottleneck/momentsld_${MODE}.json" \
  --restarts "${RESTARTS:-32}" \
  --seed "${SEED:-7411}" \
  "${EXTRA[@]}"
