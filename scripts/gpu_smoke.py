from __future__ import annotations

import json

import torch

from ghist_bottleneck.model import BottleneckNet, mixture_nll


def main() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable in the scheduled GPU job")
    device = torch.device("cuda")
    model = BottleneckNet().to(device)
    windows = torch.randn(16, 64, 8, device=device)
    context = torch.randn(16, 116, device=device)
    target = torch.randn(16, 2, device=device)
    logits, means, scale = model.mixture_parameters(windows, context)
    loss = mixture_nll(logits, means, scale, target)
    loss.backward()
    torch.cuda.synchronize()
    report = {
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "device": torch.cuda.get_device_name(0),
        "loss": float(loss.detach().cpu()),
        "peak_memory_bytes": int(torch.cuda.max_memory_allocated()),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

