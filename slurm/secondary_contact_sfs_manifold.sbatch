#!/usr/bin/env bash
#SBATCH --job-name=sc-sfs-manifold
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
.venv/bin/python scripts/diagnose_secondary_contact_sfs_manifold.py \
  --vcf data/raw/demography/secondary_contact/testing/secondary_contact.testing.vcf.gz \
  --chunks 'data/processed/secondary_contact/generator_sfs_chunks/chunk_*.npz' \
  --out artifacts/task2_secondary_contact/generator_sfs_manifold_diagnostic.json \
  --reference 32400 19765 68000 4770 1748 0.0001982
