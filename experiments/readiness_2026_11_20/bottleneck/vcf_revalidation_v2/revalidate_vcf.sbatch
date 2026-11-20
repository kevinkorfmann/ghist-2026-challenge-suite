#!/usr/bin/env bash
#SBATCH --job-name=ghist-task1-vcf-check
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-bottleneck-v2-20261120/vcf-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/readiness-bottleneck-v2-20261120"
export PYTHONPATH="$PWD/src:$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="$PWD/.mplconfig"
OLD="$ROOT/reviews/readiness-bottleneck-20261120"
"$ROOT/.venv/bin/python" -m scripts.validate_bottleneck_readiness \
  --out "$PWD/results" --case 0 --vcf-only
"$ROOT/.venv/bin/python" -m scripts.revalidate_bottleneck_vcf \
  --trees "$OLD/case12-revalidation-verified.trees" \
  --spectrum "$OLD/results/case_12/spectrum.npz" --protocol "$OLD/results/protocol.json" \
  --out "$PWD/case12_vcf_contract.json"
