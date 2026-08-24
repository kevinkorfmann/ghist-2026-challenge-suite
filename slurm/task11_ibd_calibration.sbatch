#!/usr/bin/env bash
#SBATCH --job-name=task11-ibd-calibration
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
PED_ROOT="$PROJECT_ROOT/external/ped-sim"
MAP_ROOT="$PROJECT_ROOT/external/dog-recombination"
OUT_ROOT="$PROJECT_ROOT/artifacts/task11_multifamily_inbreeding/calibration"
mkdir -p "$OUT_ROOT"

if [[ ! -s "$MAP_ROOT/campbell2016_average.simmap" ]]; then
  cd "$PROJECT_ROOT"
  PYTHONPATH=src .venv/bin/python scripts/plink_map_to_simmap.py \
    --input "$MAP_ROOT/campbell2016_average.plink.map" \
    --output "$MAP_ROOT/campbell2016_average.simmap"
fi

seed=$((91126000 + SLURM_ARRAY_TASK_ID))
"$PED_ROOT/ped-sim" \
  -d "$PROJECT_ROOT/configs/task11_pedigrees.def" \
  -m "$MAP_ROOT/campbell2016_average.simmap" \
  -o "$OUT_ROOT/replicate_${SLURM_ARRAY_TASK_ID}" \
  --intf "$PED_ROOT/interfere/nu_p_campbell.tsv" \
  --seed "$seed"
