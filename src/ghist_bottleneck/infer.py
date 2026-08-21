from __future__ import annotations

from pathlib import Path

import numpy as np
import torch
import yaml
import joblib

from .features import FeatureBundle
from .model import BottleneckNet
from .train import reference_features


def predict(observed: str | Path, checkpoint_dir: str | Path, out: str | Path) -> dict:
    b = FeatureBundle.load(observed)
    checkpoint_dir = Path(checkpoint_dir)
    norm = np.load(checkpoint_dir / "normalization.npz")
    xw = ((b.windows[None] - norm["mean"]) / norm["std"]).astype(np.float32)
    xs = b.sfs[None].astype(np.float32)
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    preds, aleatoric = [], []
    for path in sorted(checkpoint_dir.glob("model_*.pt")):
        model = BottleneckNet(b.windows.shape[1], b.sfs.shape[0]).to(device)
        state = torch.load(path, map_location=device)
        model.load_state_dict(state["state_dict"])
        model.eval()
        with torch.no_grad():
            mean, log_var = model(torch.from_numpy(xw).to(device), torch.from_numpy(xs).to(device))
        preds.append(mean.cpu().numpy()[0])
        aleatoric.append(np.exp(log_var.cpu().numpy()[0]))
    if not preds:
        raise FileNotFoundError(f"No model_*.pt checkpoints in {checkpoint_dir}")
    z = np.asarray(preds)
    deep_mean = z.mean(0) * norm["target_std"][0] + norm["target_mean"][0]
    total_var = (np.asarray(aleatoric).mean(0) + z.var(0)) * np.square(norm["target_std"][0])
    mean = deep_mean
    component_predictions = {"deep_ensemble_latent": deep_mean.tolist()}
    reference_path = checkpoint_dir / "reference.joblib"
    if reference_path.exists():
        reference = joblib.load(reference_path)
        raw_w = b.windows[None].astype(np.float32)
        ref_mean = reference["model"].predict(reference_features(raw_w, xs))[0]
        weights = np.asarray(reference["weights"])
        mean = weights[0] * deep_mean + weights[1] * ref_mean
        component_predictions["reference_latent"] = ref_mean.tolist()
        component_predictions["blend_weights_deep_reference"] = weights.T.tolist()
    generations = float(np.exp(mean[0]))
    fraction = float(1 / (1 + np.exp(-mean[1])))
    result = {
        "generations": int(round(generations)),
        "post_decline_fraction": round(fraction, 6),
        "diagnostics": {
            "ensemble_members": len(preds),
            "latent_sd": [float(x) for x in np.sqrt(total_var)],
            "n_variants": b.n_variants,
            "n_haplotypes": b.n_haplotypes,
            "components": component_predictions,
        },
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump({k: result[k] for k in ("generations", "post_decline_fraction")}, sort_keys=False))
    return result
