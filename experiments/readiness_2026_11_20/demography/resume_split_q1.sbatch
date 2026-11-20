#!/bin/bash
#SBATCH --job-name=ghist-resume-split
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/2026-11-20/artifacts/readiness/demography/split-q1-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$ROOT/reviews/2026-11-20"
export PYTHONPATH="$PWD/src:$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/.mplconfig" XDG_CACHE_HOME="$ROOT/.cache"
"$ROOT/.venv/bin/python" -m scripts.simulate_split_migration_bgs_selection \
 --coding-bed "$ROOT/data/raw/shared/demography/BGS_demography_challenge_coding_regions.bed" \
 --out "$PWD/artifacts/readiness/demography/split_q1.json" \
 --trees "$PWD/artifacts/scientific/demography/split_q1.trees" \
 --resume --existing-trees-sha256 30c740212a73e15cf8e36eef44470bfa8890198e6192f7805791f9db6f447eb4 \
 --slim-path "$ROOT/external/SLiM-v4.3/build/slim" --ancestral-size 11000 \
 --afr-size 8000 --eur-size 19000 --generations 7800 --migration-rate 0.00007 \
 --left 5000000 --right 6000000 --afr-samples 12 --eur-samples 10 \
 --slim-scaling-factor 1 --slim-burn-in 2 --fit-restarts 12 --seed 202611221
