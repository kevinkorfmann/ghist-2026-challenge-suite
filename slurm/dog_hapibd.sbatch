#!/usr/bin/env bash
#SBATCH --job-name=dog-hapibd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4096M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TOOL_ROOT="$PROJECT_ROOT/external/hap-ibd"
MAP_ROOT="$PROJECT_ROOT/external/dog-recombination"
OUT_ROOT="$PROJECT_ROOT/artifacts/task11_multifamily_inbreeding/hapibd"
mkdir -p "$TOOL_ROOT" "$MAP_ROOT" "$OUT_ROOT"

if [[ ! -s "$TOOL_ROOT/hap-ibd.jar" ]]; then
  wget -q https://faculty.washington.edu/browning/hap-ibd.jar -O "$TOOL_ROOT/hap-ibd.jar"
fi
if [[ ! -s "$MAP_ROOT/dog_genetic_maps.tar.gz" ]]; then
  wget -q \
    https://raw.githubusercontent.com/cflerin/dog_recombination/master/dog_genetic_maps.tar.gz \
    -O "$MAP_ROOT/dog_genetic_maps.tar.gz"
fi
if [[ ! -f "$MAP_ROOT/chr1_average_canFam3.1.txt" ]]; then
  tar -xzf "$MAP_ROOT/dog_genetic_maps.tar.gz" -C "$MAP_ROOT"
fi
cd "$PROJECT_ROOT"
PYTHONPATH=src .venv/bin/python scripts/dog_maps_to_plink.py \
  --map-dir "$MAP_ROOT" \
  --output "$MAP_ROOT/campbell2016_average.plink.map"

PREPARED="$OUT_ROOT/input.major_imputed.vcf.gz"
PYTHONPATH=src .venv/bin/python scripts/prepare_hapibd_vcf.py \
  --input data/raw/relatedness/multifamily_inbreeding/testing/relatedness3.testing.vcf.gz \
  --output "$PREPARED"

java -Xmx14g -jar "$TOOL_ROOT/hap-ibd.jar" \
  gt="$PREPARED" \
  map="$MAP_ROOT/campbell2016_average.plink.map" \
  out="$OUT_ROOT/segments" \
  min-seed=1.0 \
  min-extend=0.2 \
  min-output=2.0 \
  min-markers=50 \
  max-gap=1000 \
  nthreads="$SLURM_CPUS_PER_TASK"
