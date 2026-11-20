#!/usr/bin/env bash
#SBATCH --job-name=human-hapibd
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
TOOL_ROOT="$PROJECT_ROOT/external/hap-ibd"
MAP_ROOT="$PROJECT_ROOT/external/ped-sim-map"
mkdir -p "$TOOL_ROOT" "$PROJECT_ROOT/artifacts"

if [[ ! -s "$TOOL_ROOT/hap-ibd.jar" ]]; then
  wget -q https://faculty.washington.edu/browning/hap-ibd.jar -O "$TOOL_ROOT/hap-ibd.jar"
fi
if [[ ! -s "$MAP_ROOT/refined_average.plink.map" ]]; then
  cd "$PROJECT_ROOT"
  PYTHONPATH=src .venv/bin/python scripts/sex_map_to_plink.py \
    --input "$MAP_ROOT/refined_mf.simmap" \
    --output "$MAP_ROOT/refined_average.plink.map"
fi

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  INPUT="$PROJECT_ROOT/data/raw/relatedness/large_family/testing/relatedness2_testing.vcf.gz"
  OUT_ROOT="$PROJECT_ROOT/artifacts/task10_large_family/hapibd"
else
  INPUT="$PROJECT_ROOT/data/raw/relatedness/large_family_background_relatedness/testing/relatedness4.testing.vcf.gz"
  OUT_ROOT="$PROJECT_ROOT/artifacts/task12_large_family_background/hapibd"
fi
mkdir -p "$OUT_ROOT"
PREPARED="$OUT_ROOT/input.major_imputed.vcf.gz"
cd "$PROJECT_ROOT"
PYTHONPATH=src .venv/bin/python scripts/prepare_hapibd_vcf.py \
  --input "$INPUT" \
  --output "$PREPARED"

java -Xmx14g -jar "$TOOL_ROOT/hap-ibd.jar" \
  gt="$PREPARED" \
  map="$MAP_ROOT/refined_average.plink.map" \
  out="$OUT_ROOT/segments" \
  min-seed=1.0 \
  min-extend=0.2 \
  min-output=2.0 \
  min-markers=100 \
  max-gap=1000 \
  nthreads="$SLURM_CPUS_PER_TASK"
