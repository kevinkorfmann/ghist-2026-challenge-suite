#!/usr/bin/env bash
# Root-owned retry wrapper; this file does not submit or select a final candidate.
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:45:00
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
cd "$TASK_ROOT/source_v5"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PROJECT_ROOT/.venv/bin/python" -I - "$TASK_ROOT/stage_task11_map_bounds_v1.py" 8ef8ec5403285013a9414e4a8bfa14d12c45d109eacc7ed952b4d6996f06cd3b "$TASK_ROOT/prepare_pragmatic_context_v1.py" 23e66a96fef095709f28da9273ec5c62f2eb70f4fd00d5df163ad1fcdc11d867 <<'PYVERIFY'
import hashlib
import pathlib
import sys
for path, expected in zip(sys.argv[1::2], sys.argv[2::2], strict=True):
    if hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest() != expected:
        raise ValueError("Staging/adapter source hash mismatch")
PYVERIFY
mkdir -- "$2"
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/stage_task11_map_bounds_v1.py"   --vcf "$1" --map "$PROJECT_ROOT/external/dog-recombination/campbell2016_average.plink.map"   --map-sha256 c10ba1cb15e0f50888ad37b87fceea7dfb232b74634d08762d88bd7005cd851e   --min-markers 50 --out "$2/map_stage"
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/prepare_pragmatic_context_v1.py" "$TASK_ROOT"   --vcf "$2/map_stage/map_bounded.vcf.gz" --out "$2/inference" --dataset final   --candidate-sha256 6ccfda5f8133b5dee5265f35ac3a328fc4c35cfa88af076e11b9873af781ea15   --runner-sha256 7e2bdb9ffcd8f42802e2470a93f2c325b6857d57a963854b74b37a41ef1f04ce
