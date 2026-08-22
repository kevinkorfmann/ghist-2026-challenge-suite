#!/usr/bin/env bash
#SBATCH --job-name=fastcxt-bootstrap
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11 cuda/12.8.1

PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
FASTCXT_ROOT="$PROJECT_ROOT/external/fastcxt"
FASTCXT_REVISION="7fb89f451cca734e03bb974bca7faffe9a94e194"

mkdir -p "$PROJECT_ROOT/external" "$PROJECT_ROOT/checkpoints/fastcxt"
if [[ ! -d "$FASTCXT_ROOT/.git" ]]; then
  git clone --filter=blob:none https://github.com/kevinkorfmann/fastcxt.git "$FASTCXT_ROOT"
fi
git -C "$FASTCXT_ROOT" fetch origin "$FASTCXT_REVISION"
git -C "$FASTCXT_ROOT" checkout --detach "$FASTCXT_REVISION"

"$PROJECT_ROOT/.venv/bin/pip" install packaging ninja einops lightning
"$PROJECT_ROOT/.venv/bin/pip" install --no-deps --no-build-isolation \
  "causal-conv1d==1.6.2.post1"
"$PROJECT_ROOT/.venv/bin/pip" install --no-deps --no-build-isolation \
  "mamba-ssm==2.3.2.post1"
"$PROJECT_ROOT/.venv/bin/pip" install --no-deps -e "$FASTCXT_ROOT"

"$PROJECT_ROOT/.venv/bin/python" - <<'PY'
import importlib.metadata
import json

import torch

from fastcxt.config import FastCxtConfig
from fastcxt.model import FastCxtModel

if not torch.cuda.is_available():
    raise RuntimeError("CUDA is required for the FastCXT runtime gate")

config = FastCxtConfig(
    d_model=64,
    d_state=16,
    n_enc_layers=1,
    n_dec_layers=1,
    max_samples=16,
    n_windows=32,
    stem_channels=16,
    conv_channels=16,
    kernel_sizes=(3, 11),
    dropout=0.0,
    device="cuda",
    batch_size=2,
)
model = FastCxtModel(config).cuda().train()
source = torch.randn(2, 2, config.n_windows, 12, device="cuda")
mutation_rate = torch.full((2, 1), -17.2, device="cuda")
with torch.autocast("cuda", dtype=torch.bfloat16):
    output = model(source, mutation_rate)
    loss = output.square().mean()
loss.backward()
torch.cuda.synchronize()

print(json.dumps({
    "revision": "7fb89f451cca734e03bb974bca7faffe9a94e194",
    "torch": torch.__version__,
    "causal_conv1d": importlib.metadata.version("causal-conv1d"),
    "mamba_ssm": importlib.metadata.version("mamba-ssm"),
    "gpu": torch.cuda.get_device_name(0),
    "parameters": sum(parameter.numel() for parameter in model.parameters()),
    "output_shape": list(output.shape),
    "loss": float(loss.detach().cpu()),
}))
PY
