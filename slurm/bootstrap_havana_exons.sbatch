#!/usr/bin/env bash
#SBATCH --job-name=havana-exons
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="$PROJECT_ROOT/data/reference/sweeps"
GTF="$REFERENCE_ROOT/Homo_sapiens.GRCh38.104.chr.gtf.gz"
mkdir -p "$REFERENCE_ROOT"

if [[ ! -s "$GTF" ]]; then
  wget -q \
    https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz \
    -O "$GTF"
fi

cd "$PROJECT_ROOT"
PYTHONPATH=src .venv/bin/python scripts/extract_havana_exons.py \
  --gtf "$GTF" \
  --chromosome 7 \
  --bed data/reference/sweeps/havana104_chr7_exons.bed \
  --summary artifacts/shared_sweeps/havana104_chr7_exons_summary.json
