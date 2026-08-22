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
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


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
    model_config: dict | None = None,
    learning_rate: float = 2e-4,
    weight_decay: float = 2e-3,
    batch_size: int = 64,
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
    model_config = dict(model_config or {})
    resolved_model_config = {
        "window_features": int(xw.shape[2]),
        "sfs_bins": int(xs.shape[1]),
        **model_config,
    }
    epochs = 4 if quick else epochs
    for seed in seeds[:1] if quick else seeds:
        _seed_everything(seed)
        model = BottleneckNet(**resolved_model_config).to(device)
        opt = torch.optim.AdamW(
            model.parameters(), lr=learning_rate, weight_decay=weight_decay
        )
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
        train_ds = TensorDataset(*[torch.from_numpy(a[train_idx]).float() for a in (xw, xs, y)])
        loader = DataLoader(
            train_ds, batch_size=min(batch_size, len(train_ds)), shuffle=True
        )
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
                torch.save(
                    {
                        "state_dict": model.state_dict(),
                        "seed": seed,
                        "model_config": resolved_model_config,
                    },
                    model_path,
                )
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

    deep_preds, deep_vars = [], []
    vw = torch.from_numpy(xw[val_idx]).float().to(device)
    vs = torch.from_numpy(xs[val_idx]).float().to(device)
    for model_path in trained_paths:
        checkpoint = torch.load(model_path, map_location=device)
        model = BottleneckNet(**checkpoint["model_config"]).to(device)
        model.load_state_dict(checkpoint["state_dict"])
        model.eval()
        with torch.no_grad():
            pred, log_var = model(vw, vs)
        pred = pred.cpu().numpy() * target_std + target_mean
        variance = np.exp(log_var.cpu().numpy()) * np.square(target_std)
        deep_preds.append(pred)
        deep_vars.append(variance)
    deep_pred = np.mean(deep_preds, axis=0)
    total_variance = np.mean(deep_vars, axis=0) + np.var(deep_preds, axis=0)
    lower = deep_pred - 1.6448536269514722 * np.sqrt(total_variance)
    upper = deep_pred + 1.6448536269514722 * np.sqrt(total_variance)
    coverage_90 = np.mean(
        (y_raw[val_idx] >= lower) & (y_raw[val_idx] <= upper), axis=0
    )
    deep_mse = np.mean(np.square(deep_pred - y_raw[val_idx]), axis=0)
    ref_mse = np.mean(np.square(ref_pred - y_raw[val_idx]), axis=0)
    inv = np.stack((1 / (deep_mse + 1e-6), 1 / (ref_mse + 1e-6)), axis=0)
    weights = inv / inv.sum(axis=0, keepdims=True)
    scenario_report = {}
    validation_scenarios = np.zeros(len(val_idx), dtype=np.int16)
    if "nuisance" in d:
        scenarios = d["nuisance"][:, 0].astype(int)
        validation_scenarios = scenarios[val_idx].astype(np.int16)
        for scenario in np.unique(scenarios[val_idx]):
            mask = scenarios[val_idx] == scenario
            scenario_report[str(scenario)] = {
                "examples": int(mask.sum()),
                "deep_rmse_latent": np.sqrt(
                    np.mean(np.square(deep_pred[mask] - y_raw[val_idx][mask]), axis=0)
                ).tolist(),
                "reference_rmse_latent": np.sqrt(
                    np.mean(np.square(ref_pred[mask] - y_raw[val_idx][mask]), axis=0)
                ).tolist(),
                "deep_coverage_90": np.mean(
                    (y_raw[val_idx][mask] >= lower[mask])
                    & (y_raw[val_idx][mask] <= upper[mask]),
                    axis=0,
                ).tolist(),
            }
    np.savez_compressed(
        out_dir / "validation_predictions.npz",
        target_latent=y_raw[val_idx],
        deep_members_latent=np.asarray(deep_preds),
        deep_mean_latent=deep_pred,
        deep_total_variance_latent=total_variance,
        reference_latent=ref_pred,
        nuisance_scenario=validation_scenarios,
        validation_indices=val_idx,
    )
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
        "heldout_deep_coverage_90": coverage_90.tolist(),
        "validation_by_nuisance_scenario": scenario_report,
        "configuration": {
            "model": resolved_model_config,
            "training": {
                "epochs": epochs,
                "learning_rate": learning_rate,
                "weight_decay": weight_decay,
                "batch_size": batch_size,
                "seeds": list(seeds[:1] if quick else seeds),
            },
        },
    }
    (out_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
