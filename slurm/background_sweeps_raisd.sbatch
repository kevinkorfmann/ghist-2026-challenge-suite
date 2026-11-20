#!/usr/bin/env bash
#SBATCH --job-name=bgs-sweeps-raisd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="$PROJECT_ROOT/external/raisd"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  VCF="$PROJECT_ROOT/data/raw/sweeps/single_sweep_background_selection/testing/onesweep_bgs.testing.vcf.gz"
  RESULT_ROOT="$PROJECT_ROOT/artifacts/task7_single_sweep_bgs/raisd"
  LABEL="task7_w50"
  SUMMARY="artifacts/task7_single_sweep_bgs/raisd_summary.json"
else
  VCF="$PROJECT_ROOT/data/raw/sweeps/multiple_sweeps_background_selection/testing/multisweep_bgs.testing.vcf.gz"
  RESULT_ROOT="$PROJECT_ROOT/artifacts/task8_multiple_sweeps_bgs/raisd"
  LABEL="task8_w50"
  SUMMARY="artifacts/task8_multiple_sweeps_bgs/raisd_summary.json"
fi
mkdir -p "$RESULT_ROOT"

export LD_LIBRARY_PATH="$REFERENCE_ROOT/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
TASK_SCRATCH="$(mktemp -d /tmp/bgs-sweeps-raisd.XXXXXX)"
trap 'rm -rf -- "$TASK_SCRATCH"' EXIT
gzip -dc "$VCF" > "$TASK_SCRATCH/input.vcf"
cd "$RESULT_ROOT"
"$REFERENCE_ROOT/RAiSD" -n "$LABEL" -I "$TASK_SCRATCH/input.vcf" -w 50 -R -D -f

cd "$PROJECT_ROOT"
PYTHONPATH=src .venv/bin/python scripts/summarize_single_sweep_raisd.py \
  --report "$RESULT_ROOT/RAiSD_Report.$LABEL" \
  --info "$RESULT_ROOT/RAiSD_Info.$LABEL" \
  --site-report "$RESULT_ROOT/RAiSD_SiteReport.$LABEL" \
  --job-id "$SLURM_JOB_ID" \
  --out "$SUMMARY"
