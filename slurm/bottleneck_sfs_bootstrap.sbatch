#!/usr/bin/env bash
#SBATCH --job-name=bn-sfs-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
MODE="${MODE:-unfolded}"
FOLD_FLAG=()
if [[ "$MODE" == "folded" ]]; then
  FOLD_FLAG=(--folded)
elif [[ "$MODE" != "unfolded" ]]; then
  echo "MODE must be folded or unfolded" >&2
  exit 2
fi
.venv/bin/python scripts/bootstrap_bottleneck_sfs.py \
  --vcf data/processed/task1_bottleneck/bottleneck.testing.normalized.vcf.gz \
  --fit "artifacts/task1_bottleneck/moments_${MODE}.json" \
  --out "artifacts/task1_bottleneck/moments_${MODE}_block_bootstrap.json" \
  --replicates 200 \
  --workers 16 \
  "${FOLD_FLAG[@]}"
