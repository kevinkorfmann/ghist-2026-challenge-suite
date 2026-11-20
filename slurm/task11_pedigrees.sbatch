#!/usr/bin/env bash
#SBATCH --job-name=task11-pedigrees
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TOOL_ROOT="$PROJECT_ROOT/external/bonsaitree"
cd "$PROJECT_ROOT"

if [[ ! -s "$TOOL_ROOT/setup.py" ]]; then
  git clone --depth 1 https://github.com/23andMe/bonsaitree.git "$TOOL_ROOT"
fi
.venv/bin/python -m pip install -q -e "$TOOL_ROOT"

PYTHONPATH=src:scripts .venv/bin/python scripts/cluster_task11_families.py \
  --pairs artifacts/task11_multifamily_inbreeding/hapibd/summary/pair_features.csv \
  --out artifacts/task11_multifamily_inbreeding/family_clusters/summary.json

PYTHONPATH=src:scripts .venv/bin/python scripts/reconstruct_task11_pedigrees.py \
  --ibd artifacts/task11_multifamily_inbreeding/hapibd/segments.ibd.gz \
  --map external/dog-recombination/campbell2016_average.plink.map \
  --clusters artifacts/task11_multifamily_inbreeding/family_clusters/summary.json \
  --genome-cm 1978.006833142 \
  --out artifacts/task11_multifamily_inbreeding/bonsai/summary.json
