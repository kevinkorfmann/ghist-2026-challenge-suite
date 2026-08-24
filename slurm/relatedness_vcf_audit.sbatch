#!/usr/bin/env bash
#SBATCH --job-name=relatedness-vcf-audit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

case "$SLURM_ARRAY_TASK_ID" in
  0)
    VCF="data/raw/relatedness/large_family/testing/relatedness2_testing.vcf.gz"
    OUT="artifacts/task10_large_family/audit_summary.json"
    ;;
  1)
    VCF="data/raw/relatedness/multifamily_inbreeding/testing/relatedness3.testing.vcf.gz"
    OUT="artifacts/task11_multifamily_inbreeding/audit_summary.json"
    ;;
  2)
    VCF="data/raw/relatedness/large_family_background_relatedness/testing/relatedness4.testing.vcf.gz"
    OUT="artifacts/task12_large_family_background/audit_summary.json"
    ;;
esac

PYTHONPATH=src .venv/bin/python scripts/audit_relatedness_vcf.py --vcf "$VCF" --out "$OUT"
