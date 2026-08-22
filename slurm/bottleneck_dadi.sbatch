#!/usr/bin/env bash
#SBATCH --job-name=bn-dadi
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

MODE="${MODE:-unfolded}"
EXTRA=()
if [[ "$MODE" == "folded" ]]; then
  EXTRA+=(--folded)
elif [[ "$MODE" != "unfolded" ]]; then
  echo "MODE must be folded or unfolded" >&2
  exit 2
fi

.venv/bin/python -m ghist_bottleneck.dadi_reference \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --out "artifacts/task1_bottleneck/dadi_${MODE}.json" \
  --restarts "${RESTARTS:-64}" \
  --seed "${SEED:-8241}" \
  "${EXTRA[@]}"
