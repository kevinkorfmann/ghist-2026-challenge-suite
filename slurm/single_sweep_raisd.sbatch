#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-raisd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
REFERENCE_ROOT="$PROJECT_ROOT/external/raisd"
RESULT_ROOT="$PROJECT_ROOT/artifacts/task5_single_sweep/raisd"
mkdir -p "$PROJECT_ROOT/external" "$RESULT_ROOT"

if [[ ! -d "$REFERENCE_ROOT/.git" ]]; then
  git clone --depth 1 https://github.com/alachins/raisd.git "$REFERENCE_ROOT"
fi
cd "$REFERENCE_ROOT"
if [[ ! -x RAiSD ]]; then
  ./install-RAiSD-GSL.sh
fi

export LD_LIBRARY_PATH="$REFERENCE_ROOT/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
gzip -dc \
  "$PROJECT_ROOT/data/raw/sweeps/single_sweep/testing/onesweep.testing.vcf.gz" \
  > "$SLURM_TMPDIR/onesweep.testing.vcf"
cd "$RESULT_ROOT"
"$REFERENCE_ROOT/RAiSD" \
  -n task5_w50 \
  -I "$SLURM_TMPDIR/onesweep.testing.vcf" \
  -w 50 \
  -R \
  -D \
  -f
