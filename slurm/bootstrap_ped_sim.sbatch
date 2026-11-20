#!/usr/bin/env bash
#SBATCH --job-name=task9-ped-sim-setup
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load gcc/13.3.0
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
PED_ROOT="$PROJECT_ROOT/external/ped-sim"
MAP_ROOT="$PROJECT_ROOT/external/ped-sim-map"
GSL_ROOT="$PROJECT_ROOT/external/raisd/gsl"
mkdir -p "$PROJECT_ROOT/external" "$MAP_ROOT"

if [[ ! -d "$PED_ROOT/.git" ]]; then
  git clone --depth 1 https://github.com/williamslab/ped-sim.git "$PED_ROOT"
fi
test -f "$GSL_ROOT/lib/libgsl.so"
cp "$PED_ROOT/Makefile-gsl" "$PED_ROOT/Makefile"
make -C "$PED_ROOT" clean
make -C "$PED_ROOT" -j "$SLURM_CPUS_PER_TASK" \
  CPPFLAGS="-std=c++11 -Wall -DUSEGSL -O2 -I$GSL_ROOT/include" \
  CFLAGS="-Wall -DUSEGSL -O2 -I$GSL_ROOT/include" \
  LIBS="-L$GSL_ROOT/lib -Wl,-rpath,$GSL_ROOT/lib -lz -lgsl -lgslcblas"

if [[ ! -f "$MAP_ROOT/refined_mf.simmap" ]]; then
  wget -q \
    https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/raw/master/Refined_genetic_map_b37.tar.gz \
    -O "$MAP_ROOT/Refined_genetic_map_b37.tar.gz"
  tar -xzf "$MAP_ROOT/Refined_genetic_map_b37.tar.gz" -C "$MAP_ROOT"
  printf '#chr\tpos\tmale_cM\tfemale_cM\n' > "$MAP_ROOT/refined_mf.simmap"
  for chromosome in {1..22}; do
    paste \
      "$MAP_ROOT/Refined_genetic_map_b37/male_chr${chromosome}.txt" \
      "$MAP_ROOT/Refined_genetic_map_b37/female_chr${chromosome}.txt" \
      | awk -v OFS='\t' 'NR > 1 && $2 == $6 {print $1,$2,$4,$8}' \
      | sed 's/^chr//' >> "$MAP_ROOT/refined_mf.simmap"
  done
fi

"$PED_ROOT/ped-sim" --help >/dev/null 2>&1 || [[ "$?" == "1" ]]
test -s "$MAP_ROOT/refined_mf.simmap"
