#!/usr/bin/env bash
#SBATCH --job-name=growth-select-stack
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
module load gcc/13.3.0 cmake/3.31.6 python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
EXTERNAL_ROOT="$PROJECT_ROOT/external"
SLIM_ROOT="$EXTERNAL_ROOT/SLiM-v4.3"
mkdir -p "$EXTERNAL_ROOT"

if [[ ! -d "$SLIM_ROOT/.git" ]]; then
  git clone --depth 1 --branch v4.3 https://github.com/MesserLab/SLiM.git "$SLIM_ROOT"
fi
cmake -S "$SLIM_ROOT" -B "$SLIM_ROOT/build" -DBUILD_SLIMGUI=OFF
cmake --build "$SLIM_ROOT/build" --target slim -j "$SLURM_CPUS_PER_TASK"

"$PROJECT_ROOT/.venv/bin/pip" install 'stdpopsim==0.3.0' 'pyslim>=1.0.4,<1.1'
"$SLIM_ROOT/build/slim" -v
"$PROJECT_ROOT/.venv/bin/python" - <<'PY'
import pyslim
import stdpopsim

species = stdpopsim.get_species("DroMel")
print("stdpopsim", stdpopsim.__version__)
print("pyslim", pyslim.__version__)
print("map", species.get_genetic_map("ComeronCrossoverV2_dm6").id)
print("dfe", species.get_dfe("Gamma_H17").id)
PY
