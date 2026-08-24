#!/usr/bin/env bash
#SBATCH --job-name=bonsai-pedigree
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=08:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TOOL_ROOT="$PROJECT_ROOT/external/bonsaitree"
cd "$PROJECT_ROOT"

if [[ ! -s "$TOOL_ROOT/setup.py" ]]; then
  git clone --depth 1 https://github.com/23andMe/bonsaitree.git "$TOOL_ROOT"
fi
.venv/bin/python -m pip install -q -e "$TOOL_ROOT"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  ROOT="artifacts/task10_large_family/hapibd"
  THRESHOLD="10"
else
  ROOT="artifacts/task12_large_family_background/hapibd"
  THRESHOLD="20"
fi

PYTHONPATH=src .venv/bin/python scripts/reconstruct_bonsai_pedigree.py \
  --ibd "$ROOT/segments.ibd.gz" \
  --map external/ped-sim-map/refined_average.plink.map \
  --threshold-cm "$THRESHOLD" \
  --out "$ROOT/bonsai/summary.json"
