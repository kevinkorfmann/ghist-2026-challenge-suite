from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

from ghist_secondary_contact.features import transform_targets
from ghist_secondary_contact.model import SecondaryContactNet, mixture_nll


def _load(pattern: str) -> dict[str, np.ndarray]:
    paths = []
    for item in pattern.split(","):
        paths.extend(Path(path) for path in sorted(glob.glob(item.strip())))
    paths = list(dict.fromkeys(paths))
    if not paths:
        raise FileNotFoundError(f"no dataset chunks match {pattern!r}")
    required = (
        "cell_id",
        "point_id",
        "partition",
        "truth",
        "joint_sfs",
    )
    loaded = [np.load(path) for path in paths]
    # Spectrum-surrogate chunks use point/replicate identifiers instead of a
    # cell identifier and intentionally omit chromosome-context features.
    available = set(loaded[0].files)
    if "cell_id" not in available:
        required = tuple(key for key in required if key != "cell_id")
    result = {key: np.concatenate([chunk[key] for chunk in loaded]) for key in required}
    if "cell_id" in result and len(np.unique(result["cell_id"])) != len(result["cell_id"]):
        raise ValueError("dataset chunks contain duplicate cell identifiers")
    if "windows" in available:
        result["windows"] = np.concatenate([chunk["windows"] for chunk in loaded])
    else:
        result["windows"] = np.zeros((len(result["truth"]), 1, 1), dtype=np.float32)
    if "global_summary" in available:
        result["global_summary"] = np.concatenate(
            [chunk["global_summary"] for chunk in loaded]
        )
    else:
        counts = np.concatenate([chunk["n_variants"] for chunk in loaded])
        result["global_summary"] = np.log1p(counts.astype(np.float64))[:, None]
    result["design_domain"] = np.concatenate(
        [
            chunk["design_domain"]
            if "design_domain" in chunk.files
            else np.full(len(chunk["truth"]), 2, dtype=np.int8)
            for chunk in loaded
        ]
    )
    return result


def _point_rrmse(
    prediction: np.ndarray,
    truth: np.ndarray,
    point_ids: np.ndarray,
) -> float:
    point_prediction, point_truth = [], []
    for point_id in np.unique(point_ids):
        selected = point_ids == point_id
        point_prediction.append(prediction[selected].mean(axis=0))
        point_truth.append(truth[selected][0])
    relative = (np.asarray(point_prediction) - np.asarray(point_truth)) / np.asarray(point_truth)
    return float(np.sqrt(np.mean(np.sum(relative**2, axis=1))))


def _physical(values: torch.Tensor, mean: torch.Tensor, scale: torch.Tensor) -> torch.Tensor:
    transformed = values * scale + mean
    split = torch.exp(transformed[:, 2].clamp(3.5, 13.0))
    fraction = torch.sigmoid(transformed[:, 3].clamp(-10.0, 10.0))
    return torch.stack(
        (
            torch.exp(transformed[:, 0].clamp(6.0, 14.5)),
            torch.exp(transformed[:, 1].clamp(6.0, 14.5)),
            split,
            split * fraction,
            torch.exp(transformed[:, 4].clamp(-20.0, -2.0)),
        ),
        dim=1,
    )


@torch.no_grad()
def _evaluate(
    model: SecondaryContactNet,
    arrays: dict[str, torch.Tensor],
    indices: np.ndarray,
    target_mean: torch.Tensor,
    target_scale: torch.Tensor,
    *,
    batch_size: int,
) -> tuple[float, np.ndarray]:
    model.eval()
    predictions = []
    for start in range(0, len(indices), batch_size):
        rows = torch.as_tensor(
            indices[start : start + batch_size], dtype=torch.long, device=target_mean.device
        )
        mean, _ = model(
            arrays["windows"][rows],
            arrays["joint_sfs"][rows, None],
            arrays["global_summary"][rows],
        )
        predictions.append(_physical(mean, target_mean, target_scale).cpu().numpy())
    prediction = np.concatenate(predictions)
    all_rows = torch.as_tensor(indices, dtype=torch.long, device=target_mean.device)
    truth = arrays["truth"][all_rows].cpu().numpy()
    point_ids = arrays["point_id"][all_rows].cpu().numpy()
    return _point_rrmse(prediction, truth, point_ids), prediction


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunks", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--batch-size", type=int, default=24)
    parser.add_argument("--learning-rate", type=float, default=2e-4)
    parser.add_argument("--weight-decay", type=float, default=2e-4)
    parser.add_argument("--score-weight", type=float, default=0.15)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--width", type=int, default=128)
    parser.add_argument("--context-backend", choices=("transformer", "mamba"), default="transformer")
    parser.add_argument("--seed", type=int, default=7701)
    parser.add_argument("--domain-filter", type=str)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = _load(args.chunks)
    if args.domain_filter:
        domains = np.asarray(
            [
                int(value)
                for value in args.domain_filter.replace(":", ",").split(",")
            ],
            dtype=np.int8,
        )
        keep = np.isin(data["design_domain"], domains)
        if not np.any(keep):
            raise ValueError(f"domain filter retained no rows: {domains.tolist()}")
        data = {key: value[keep] for key, value in data.items()}
    calibration = data["partition"] == 0
    if not np.any(calibration) or not np.any(data["partition"] == 1):
        raise ValueError("calibration and selection partitions must both be present")
    window_mean = data["windows"][calibration].mean(axis=(0, 1), keepdims=True)
    window_scale = data["windows"][calibration].std(axis=(0, 1), keepdims=True).clip(1e-6)
    global_mean = data["global_summary"][calibration].mean(axis=0, keepdims=True)
    global_scale = data["global_summary"][calibration].std(axis=0, keepdims=True).clip(1e-6)
    transformed = transform_targets(data["truth"])
    target_mean_np = transformed[calibration].mean(axis=0)
    target_scale_np = transformed[calibration].std(axis=0).clip(1e-6)
    arrays = {
        "windows": torch.as_tensor(
            (data["windows"] - window_mean) / window_scale, dtype=torch.float32, device=device
        ),
        "joint_sfs": torch.as_tensor(data["joint_sfs"], dtype=torch.float32, device=device),
        "global_summary": torch.as_tensor(
            (data["global_summary"] - global_mean) / global_scale,
            dtype=torch.float32,
            device=device,
        ),
        "target": torch.as_tensor(
            (transformed - target_mean_np) / target_scale_np,
            dtype=torch.float32,
            device=device,
        ),
        "truth": torch.as_tensor(data["truth"], dtype=torch.float32, device=device),
        "point_id": torch.as_tensor(data["point_id"], dtype=torch.long, device=device),
    }
    target_mean = torch.as_tensor(target_mean_np, dtype=torch.float32, device=device)
    target_scale = torch.as_tensor(target_scale_np, dtype=torch.float32, device=device)
    train_indices = np.flatnonzero(calibration)
    selection_indices = np.flatnonzero(data["partition"] == 1)
    evaluation_indices = np.flatnonzero(data["partition"] == 2)
    loader = DataLoader(
        TensorDataset(torch.as_tensor(train_indices, dtype=torch.long)),
        batch_size=args.batch_size,
        shuffle=True,
        generator=torch.Generator().manual_seed(args.seed),
    )
    model = SecondaryContactNet(
        window_features=data["windows"].shape[-1],
        global_features=data["global_summary"].shape[-1],
        width=args.width,
        context_backend=args.context_backend,
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    best_score = np.inf
    best_epoch = 0
    history = []
    for epoch in range(args.epochs):
        model.train()
        epoch_loss = 0.0
        for (rows_cpu,) in loader:
            rows = rows_cpu.to(device)
            logits, means, scale = model.mixture_parameters(
                arrays["windows"][rows],
                arrays["joint_sfs"][rows, None],
                arrays["global_summary"][rows],
            )
            posterior_loss = mixture_nll(logits, means, scale, arrays["target"][rows])
            weights = torch.softmax(logits, dim=1)
            transformed_mean = torch.sum(weights[..., None] * means, dim=1)
            physical = _physical(transformed_mean, target_mean, target_scale)
            # Preserve the competition geometry while capping extreme random-
            # initialization ratios that would otherwise dominate posterior learning.
            relative = (physical / arrays["truth"][rows] - 1.0).clamp(-3.0, 3.0)
            relative_loss = torch.mean(torch.sum(relative.square(), dim=1))
            loss = posterior_loss + args.score_weight * relative_loss
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            epoch_loss += float(loss.detach()) * len(rows)
        scheduler.step()
        selection_score, _ = _evaluate(
            model,
            arrays,
            selection_indices,
            target_mean,
            target_scale,
            batch_size=args.batch_size,
        )
        record = {
            "epoch": epoch + 1,
            "training_loss": epoch_loss / len(train_indices),
            "selection_rrmse": selection_score,
        }
        history.append(record)
        if selection_score < best_score:
            best_score = selection_score
            best_epoch = epoch + 1
            torch.save(
                {
                    "model": model.state_dict(),
                    "model_configuration": {
                        "window_features": data["windows"].shape[-1],
                        "global_features": data["global_summary"].shape[-1],
                        "width": args.width,
                        "context_backend": args.context_backend,
                    },
                    "normalization": {
                        "window_mean": window_mean,
                        "window_scale": window_scale,
                        "global_mean": global_mean,
                        "global_scale": global_scale,
                        "target_mean": target_mean_np,
                        "target_scale": target_scale_np,
                    },
                    "epoch": epoch + 1,
                    "selection_rrmse": best_score,
                },
                args.out,
            )
        print(json.dumps(record), flush=True)
        if epoch + 1 - best_epoch >= args.patience:
            print(
                json.dumps(
                    {
                        "early_stop_epoch": epoch + 1,
                        "best_epoch": best_epoch,
                        "patience": args.patience,
                    }
                ),
                flush=True,
            )
            break
    checkpoint = torch.load(args.out, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model"])
    evaluation_score, evaluation_prediction = _evaluate(
        model,
        arrays,
        evaluation_indices,
        target_mean,
        target_scale,
        batch_size=args.batch_size,
    )
    report = {
        "best_epoch": checkpoint["epoch"],
        "selection_rrmse": checkpoint["selection_rrmse"],
        "untouched_evaluation_rrmse": evaluation_score,
        "seed": args.seed,
        "context_backend": args.context_backend,
        "domain_filter": args.domain_filter,
        "history": history,
    }
    evaluation_domains = data["design_domain"][evaluation_indices]
    evaluation_truth = data["truth"][evaluation_indices]
    evaluation_point_ids = data["point_id"][evaluation_indices]
    report["untouched_evaluation_rrmse_by_domain"] = {
        label: _point_rrmse(
            evaluation_prediction[evaluation_domains == domain],
            evaluation_truth[evaluation_domains == domain],
            evaluation_point_ids[evaluation_domains == domain],
        )
        for domain, label in ((0, "core"), (1, "guard"), (2, "broad"))
        if np.any(evaluation_domains == domain)
    }
    args.out.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".evaluation.npz"),
        indices=evaluation_indices,
        prediction=evaluation_prediction,
        truth=data["truth"][evaluation_indices],
        point_id=data["point_id"][evaluation_indices],
    )
    print(json.dumps({key: value for key, value in report.items() if key != "history"}, indent=2))


if __name__ == "__main__":
    main()
