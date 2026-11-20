#!/bin/bash
#SBATCH --job-name=ghist-growth-recovery
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/growth-recovery-%A-%a.out
#SBATCH --error=logs/growth-recovery-%A-%a.err
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/2026-11-20"
export PYTHONPATH="$PWD/src:$PWD/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export XDG_CACHE_HOME="$ROOT/.cache"
conventions=(exponential sudden)
convention=${conventions[$SLURM_ARRAY_TASK_ID]}
"$ROOT/.venv/bin/python" scripts/simulate_growth_bgs_selection.py --coding-bed "$ROOT/data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed" --out "artifacts/reliability/growth/$convention.json" --trees "artifacts/reliability/growth/$convention.trees" --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --chromosome 2L --left 5000000 --right 6000000 --ancestral-size 8000 --generations 7600 --scaling 3.2 --slim-scaling-factor 10 --slim-burn-in 1 --samples 20 --seed "$((20261120 + SLURM_ARRAY_TASK_ID))" --fit-restarts 12 --convention "$convention"
