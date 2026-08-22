#!/usr/bin/env bash
#SBATCH --job-name=reference-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"

.venv/bin/pip install "moments-popgen==1.6.0" "dadi==2.4.4"
.venv/bin/pip check
.venv/bin/python - <<'PY'
import importlib.metadata
import json

import dadi
import moments

print(json.dumps({
    "moments_popgen": importlib.metadata.version("moments-popgen"),
    "dadi": importlib.metadata.version("dadi"),
    "moments_has_two_epoch": hasattr(moments.Demographics1D, "two_epoch"),
    "dadi_has_two_epoch": hasattr(dadi.Demographics1D, "two_epoch"),
}))
PY
