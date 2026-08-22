#!/usr/bin/env bash
#SBATCH --job-name=sc-generator-recovery
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
.venv/bin/python scripts/fit_secondary_contact_spectrum_archive.py \
  --archive artifacts/task2_secondary_contact/generator_spectrum_validation.sfs.npz \
  --out "artifacts/task2_secondary_contact/generator_recovery_${CONVENTION:-generator_matched}.json" \
  --convention "${CONVENTION:-generator_matched}" \
  --projection "${PROJECTION1:-40}" "${PROJECTION2:-30}" \
  --restarts "${RESTARTS:-24}" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --maxiter "${MAXITER:-500}" \
  --seed "${SEED:-20260911}"
