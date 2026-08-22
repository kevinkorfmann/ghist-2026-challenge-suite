#!/usr/bin/env bash
#SBATCH --job-name=sc-tabular-predict
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
.venv/bin/python scripts/predict_secondary_contact_tabular.py \
  --vcf data/raw/demography/secondary_contact/testing/secondary_contact.testing.vcf.gz \
  --model "$MODEL" \
  --out "$OUT"
