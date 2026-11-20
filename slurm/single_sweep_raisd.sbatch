#!/usr/bin/env bash
#SBATCH --job-name=single-sweep-raisd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
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
  if [[ ! -f gsl/lib/libgsl.so ]]; then
    if [[ ! -f gsl-latest.tar.gz ]]; then
      wget https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    fi
    if [[ -z "$(find . -maxdepth 1 -type d -name 'gsl-*' -print -quit)" ]]; then
      tar -xzf gsl-latest.tar.gz
    fi
    GSL_SOURCE="$(find . -maxdepth 1 -type d -name 'gsl-*' | sort -V | tail -1)"
    cd "$GSL_SOURCE"
    ./configure --prefix="$REFERENCE_ROOT/gsl"
    make -j "$SLURM_CPUS_PER_TASK"
    make install
    cd "$REFERENCE_ROOT"
  fi
  cp makefiles/Makefile.GSL Makefile
  make clean || true
  make -j "$SLURM_CPUS_PER_TASK"
fi

export LD_LIBRARY_PATH="$REFERENCE_ROOT/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
TASK_SCRATCH="$(mktemp -d /tmp/single-sweep-raisd.XXXXXX)"
trap 'rm -rf -- "$TASK_SCRATCH"' EXIT
gzip -dc \
  "$PROJECT_ROOT/data/raw/sweeps/single_sweep/testing/onesweep.testing.vcf.gz" \
  > "$TASK_SCRATCH/onesweep.testing.vcf"
cd "$RESULT_ROOT"
"$REFERENCE_ROOT/RAiSD" \
  -n task5_w50 \
  -I "$TASK_SCRATCH/onesweep.testing.vcf" \
  -w 50 \
  -R \
  -D \
  -f
