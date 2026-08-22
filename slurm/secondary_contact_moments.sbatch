#!/usr/bin/env bash
#SBATCH --job-name=sc-moments
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
MODE="${MODE:-folded}"
EXTRA=()
if [[ "$MODE" != *unfolded* ]]; then EXTRA+=(--folded); fi
if [[ -n "${PROJECTION1:-}" && -n "${PROJECTION2:-}" ]]; then
  EXTRA+=(--projection "$PROJECTION1" "$PROJECTION2")
elif [[ -n "${PROJECTION:-}" ]]; then
  read -r PROJECTION1 PROJECTION2 <<< "$PROJECTION"
  EXTRA+=(--projection "$PROJECTION1" "$PROJECTION2")
fi
.venv/bin/python scripts/run_secondary_contact_moments.py \
  --vcf data/raw/demography/secondary_contact/testing/secondary_contact.testing.vcf.gz \
  --out "artifacts/task2_secondary_contact/moments_${MODE}.json" \
  --restarts "${RESTARTS:-64}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --maxiter "${MAXITER:-750}" \
  --convention "${CONVENTION:-simple}" \
  "${EXTRA[@]}"
