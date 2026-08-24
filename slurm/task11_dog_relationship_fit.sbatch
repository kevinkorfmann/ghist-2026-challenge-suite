#!/usr/bin/env bash
#SBATCH --job-name=task11-dog-fit
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

thresholds=(2 5 10)
threshold="${thresholds[$SLURM_ARRAY_TASK_ID]}"

PYTHONPATH=src:scripts .venv/bin/python scripts/fit_task11_dog_relationships.py \
  --calibration \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_0.seg \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_1.seg \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_2.seg \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_3.seg \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_4.seg \
    artifacts/task11_multifamily_inbreeding/calibration/replicate_5.seg \
  --validation artifacts/task11_multifamily_inbreeding/calibration/replicate_6.seg \
  --evaluation artifacts/task11_multifamily_inbreeding/calibration/replicate_7.seg \
  --observed-ibd artifacts/task11_multifamily_inbreeding/hapibd/segments.ibd.gz \
  --observed-hbd artifacts/task11_multifamily_inbreeding/hapibd/segments.hbd.gz \
  --hbd-threshold-cm "$threshold" \
  --map external/dog-recombination/campbell2016_average.plink.map \
  --clusters artifacts/task11_multifamily_inbreeding/family_clusters/summary.json \
  --out-dir "artifacts/task11_multifamily_inbreeding/pair_model_hbd${threshold}"
