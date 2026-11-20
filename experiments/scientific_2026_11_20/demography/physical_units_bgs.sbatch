#!/bin/bash
#SBATCH --job-name=ghist-physical-bgs
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-3
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/2026-11-20/artifacts/scientific/demography/bgs-%A_%a.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/2026-11-20"
export PYTHONPATH="$PWD/src:$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/.mplconfig" XDG_CACHE_HOME="$ROOT/.cache"
OUT="$PWD/artifacts/scientific/demography"
SCALE=5
if (( SLURM_ARRAY_TASK_ID % 2 == 1 )); then SCALE=10; fi
if (( SLURM_ARRAY_TASK_ID < 2 )); then
  "$ROOT/.venv/bin/python" -m scripts.simulate_growth_bgs_selection \
    --coding-bed "$ROOT/data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed" \
    --out "$OUT/growth_q${SCALE}.json" --trees "$OUT/growth_q${SCALE}.trees" \
    --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --ancestral-size 7000 \
    --generations 9200 --scaling 3.7 --convention exponential --left 5000000 --right 6000000 \
    --samples 20 --slim-scaling-factor "$SCALE" --slim-burn-in 2 --fit-restarts 12 --seed "$((202611210 + SLURM_ARRAY_TASK_ID))"
else
  "$ROOT/.venv/bin/python" -m scripts.simulate_split_migration_bgs_selection \
    --coding-bed "$ROOT/data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed" \
    --out "$OUT/split_q${SCALE}.json" --trees "$OUT/split_q${SCALE}.trees" \
    --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --ancestral-size 11000 \
    --afr-size 8000 --eur-size 19000 --generations 7800 --migration-rate 0.00007 \
    --left 5000000 --right 6000000 --afr-samples 12 --eur-samples 10 \
    --slim-scaling-factor "$SCALE" --slim-burn-in 2 --fit-restarts 12 --seed "$((202611210 + SLURM_ARRAY_TASK_ID))"
fi
