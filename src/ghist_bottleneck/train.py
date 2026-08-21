from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import torch
import joblib
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import GroupShuffleSplit
from torch.utils.data import DataLoader, TensorDataset

from .model import BottleneckNet, mixture_nll


def _seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def _device() -> torch.device:
    return torch.device("mps" if torch.backends.mps.is_available() else "cpu")


def reference_features(windows: np.ndarray, sfs: np.ndarray) -> np.ndarray:
    """Permutation-invariant summaries for the independent tree estimator."""
    pieces = [windows.mean(1), windows.std(1)]
    pieces.extend(np.quantile(windows, q, axis=1) for q in (0.1, 0.25, 0.5, 0.75, 0.9))
    return np.concatenate([*pieces, sfs], axis=1)


def train_ensemble(
    dataset: str | Path,
    out_dir: str | Path,
    *,
    seeds: tuple[int, ...] = (17, 29, 43, 71, 101),
    epochs: int = 80,
    quick: bool = False,
) -> None:
    d = np.load(dataset)
    xw, xs, y, ne = d["windows"], d["sfs"], d["targets"], d["ancestral_ne"]
    y_raw = y.copy()
    # Hold out joint parameter cells. This is stricter than a random simulation
    # split and exposes interpolation failures that would otherwise look good.
    g_gen = np.digitize(y[:, 0], np.quantile(y[:, 0], np.linspace(0, 1, 7)[1:-1]))
    g_frac = np.digitize(y[:, 1], np.quantile(y[:, 1], np.linspace(0, 1, 6)[1:-1]))
    g_ne = np.digitize(np.log(ne), np.quantile(np.log(ne), np.linspace(0, 1, 5)[1:-1]))
    groups = g_gen * 100 + g_frac * 10 + g_ne
    train_idx, val_idx = next(GroupShuffleSplit(test_size=0.2, random_state=991).split(xw, y, groups))
    mu = xw[train_idx].mean((0, 1), keepdims=True)
    sigma = xw[train_idx].std((0, 1), keepdims=True) + 1e-6
    xw = (xw - mu) / sigma
    target_mean = y[train_idx].mean(0, keepdims=True)
    target_std = y[train_idx].std(0, keepdims=True) + 1e-6
    y = (y - target_mean) / target_std
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_dir / "normalization.npz",
        mean=mu,
        std=sigma,
        target_mean=target_mean,
        target_std=target_std,
    )
    device = _device()
    metrics = []
    trained_paths = []
    epochs = 4 if quick else epochs
    for seed in seeds[:1] if quick else seeds:
        _seed_everything(seed)
        model = BottleneckNet(xw.shape[2], xs.shape[1]).to(device)
        opt = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=2e-3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
        train_ds = TensorDataset(*[torch.from_numpy(a[train_idx]).float() for a in (xw, xs, y)])
        loader = DataLoader(train_ds, batch_size=min(64, len(train_ds)), shuffle=True)
        best = float("inf")
        for _ in range(epochs):
            model.train()
            for bw, bs, by in loader:
                bw, bs, by = bw.to(device), bs.to(device), by.to(device)
                opt.zero_grad(set_to_none=True)
                logits, means, scale = model.mixture_parameters(bw, bs)
                loss = mixture_nll(logits, means, scale, by)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                opt.step()
            scheduler.step()
            model.eval()
            with torch.no_grad():
                vw = torch.from_numpy(xw[val_idx]).float().to(device)
                vs = torch.from_numpy(xs[val_idx]).float().to(device)
                vy = torch.from_numpy(y[val_idx]).float().to(device)
                logits, means, scale = model.mixture_parameters(vw, vs)
                val = float(mixture_nll(logits, means, scale, vy).cpu())
            if val < best:
                best = val
                model_path = out_dir / f"model_{seed}.pt"
                torch.save({"state_dict": model.state_dict(), "seed": seed}, model_path)
        metrics.append({"seed": seed, "validation_nll": best})
        trained_paths.append(model_path)

    # Independent, permutation-invariant reference model. It cannot exploit the
    # ordering cues available to the conv-transformer, so agreement is meaningful.
    ref = ExtraTreesRegressor(
        n_estimators=80 if quick else 600,
        min_samples_leaf=2 if quick else 4,
        max_features=0.7,
        n_jobs=-1,
        random_state=4401,
    )
    ref_x = reference_features(d["windows"], xs)
    ref.fit(ref_x[train_idx], y_raw[train_idx])
    ref_pred = ref.predict(ref_x[val_idx])

    deep_preds = []
    vw = torch.from_numpy(xw[val_idx]).float().to(device)
    vs = torch.from_numpy(xs[val_idx]).float().to(device)
    for model_path in trained_paths:
        model = BottleneckNet(xw.shape[2], xs.shape[1]).to(device)
        model.load_state_dict(torch.load(model_path, map_location=device)["state_dict"])
        model.eval()
        with torch.no_grad():
            pred, _ = model(vw, vs)
        pred = pred.cpu().numpy() * target_std + target_mean
        deep_preds.append(pred)
    deep_pred = np.mean(deep_preds, axis=0)
    deep_mse = np.mean(np.square(deep_pred - y_raw[val_idx]), axis=0)
    ref_mse = np.mean(np.square(ref_pred - y_raw[val_idx]), axis=0)
    inv = np.stack((1 / (deep_mse + 1e-6), 1 / (ref_mse + 1e-6)), axis=0)
    weights = inv / inv.sum(axis=0, keepdims=True)
    joblib.dump(
        {"model": ref, "weights": weights, "deep_mse": deep_mse, "reference_mse": ref_mse},
        out_dir / "reference.joblib",
    )
    report = {
        "members": metrics,
        "heldout_deep_rmse_latent": np.sqrt(deep_mse).tolist(),
        "heldout_reference_rmse_latent": np.sqrt(ref_mse).tolist(),
        "blend_weights_deep_reference": weights.T.tolist(),
        "validation_examples": int(len(val_idx)),
    }
    (out_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
