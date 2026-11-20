#!/usr/bin/env bash
#SBATCH --job-name=soft-sweep-raisd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="$PROJECT_ROOT/external/raisd"
RESULT_ROOT="$PROJECT_ROOT/artifacts/task6_soft_sweep/raisd"
mkdir -p "$RESULT_ROOT"

export LD_LIBRARY_PATH="$REFERENCE_ROOT/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
TASK_SCRATCH="$(mktemp -d /tmp/soft-sweep-raisd.XXXXXX)"
trap 'rm -rf -- "$TASK_SCRATCH"' EXIT
gzip -dc \
  "$PROJECT_ROOT/data/raw/sweeps/soft_sweep/testing/softsweep.testing.vcf.gz" \
  > "$TASK_SCRATCH/softsweep.testing.vcf"
cd "$RESULT_ROOT"
"$REFERENCE_ROOT/RAiSD" \
  -n task6_w50 \
  -I "$TASK_SCRATCH/softsweep.testing.vcf" \
  -w 50 \
  -R \
  -D \
  -f

cd "$PROJECT_ROOT"
PYTHONPATH=src .venv/bin/python scripts/summarize_single_sweep_raisd.py \
  --report "$RESULT_ROOT/RAiSD_Report.task6_w50" \
  --info "$RESULT_ROOT/RAiSD_Info.task6_w50" \
  --site-report "$RESULT_ROOT/RAiSD_SiteReport.task6_w50" \
  --job-id "$SLURM_JOB_ID" \
  --out artifacts/task6_soft_sweep/raisd_summary.json
