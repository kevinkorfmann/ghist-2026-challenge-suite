#!/usr/bin/env bash
#SBATCH --job-name=cxt-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
CXT_ROOT="$PROJECT_ROOT/external/cxt"
CXT_REVISION="294b0c27226accbac830c9c922de6b2f92f83fe8"
mkdir -p "$PROJECT_ROOT/external" "$PROJECT_ROOT/checkpoints/cxt"
if [[ ! -d "$CXT_ROOT/.git" ]]; then
  git clone --filter=blob:none https://github.com/kr-colab/cxt.git "$CXT_ROOT"
fi
git -C "$CXT_ROOT" fetch origin "$CXT_REVISION"
git -C "$CXT_ROOT" checkout --detach "$CXT_REVISION"
"$PROJECT_ROOT/.venv/bin/pip" install --no-deps torchtune torchao
"$PROJECT_ROOT/.venv/bin/pip" install einops requests pandas
"$PROJECT_ROOT/.venv/bin/pip" install --no-deps -e "$CXT_ROOT"
export CXT_CHECKPOINT_CACHE="$PROJECT_ROOT/checkpoints/cxt"
"$PROJECT_ROOT/.venv/bin/python" - <<'PY'
import json
import cxt

model = cxt.load_model("w200_wmissing", device="cpu")
print(json.dumps({
    "revision": "294b0c27226accbac830c9c922de6b2f92f83fe8",
    "model": "w200_wmissing",
    "parameters": sum(parameter.numel() for parameter in model.parameters()),
}))
PY
