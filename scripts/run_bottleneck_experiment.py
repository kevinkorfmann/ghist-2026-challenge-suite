from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import yaml

from ghist_bottleneck.train import train_ensemble


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dataset_identity(path: Path) -> tuple[str, int]:
    manifest_path = path.with_suffix(path.suffix + ".manifest.json")
    size = path.stat().st_size
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if int(manifest["bytes"]) != size:
            raise ValueError("dataset size disagrees with its manifest")
        return str(manifest["sha256"]), size
    return sha256(path), size


def selection_score(output: Path, weights: dict) -> dict:
    metrics = json.loads((output / "metrics.json").read_text())
    normalization = np.load(output / "normalization.npz")
    target_scale = normalization["target_std"][0]
    overall = np.asarray(metrics["heldout_deep_rmse_latent"]) / target_scale
    nominal_metrics = metrics["validation_by_nuisance_scenario"].get("0")
    nominal = (
        np.asarray(nominal_metrics["deep_rmse_latent"]) / target_scale
        if nominal_metrics
        else overall
    )
    robust_values = [
        np.mean(np.asarray(item["deep_rmse_latent"]) / target_scale)
        for item in metrics["validation_by_nuisance_scenario"].values()
    ]
    robust_worst = max(robust_values, default=float(np.mean(overall)))
    coverage = np.asarray(metrics["heldout_deep_coverage_90"])
    calibration_gap = float(np.mean(np.abs(coverage - 0.9)))
    score = (
        float(weights["nominal_weight"]) * float(np.mean(nominal))
        + float(weights["overall_weight"]) * float(np.mean(overall))
        + float(weights["calibration_weight"]) * calibration_gap
        + float(weights["robustness_guard_weight"])
        * max(0.0, robust_worst - float(np.mean(nominal)))
    )
    return {
        "selection_score": score,
        "nominal_normalized_rmse": nominal.tolist(),
        "overall_normalized_rmse": overall.tolist(),
        "worst_scenario_mean_normalized_rmse": robust_worst,
        "coverage_90": coverage.tolist(),
        "calibration_gap": calibration_gap,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()

    registry = yaml.safe_load(args.registry.read_text())
    if args.candidate not in registry["candidates"]:
        raise KeyError(f"unknown candidate: {args.candidate}")
    configuration = registry["candidates"][args.candidate]
    output = args.out_root / args.candidate
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    training = configuration["training"]
    train_ensemble(
        args.dataset,
        output,
        model_config=configuration["model"],
        seeds=tuple(training["seeds"]),
        epochs=int(training["epochs"]),
        learning_rate=float(training["learning_rate"]),
        weight_decay=float(training["weight_decay"]),
        batch_size=int(training["batch_size"]),
    )
    revision_result = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True
    )
    if revision_result.returncode == 0:
        revision = revision_result.stdout.strip()
    else:
        marker = Path(".source-revision")
        if not marker.exists():
            raise RuntimeError("source revision is unavailable")
        revision = marker.read_text().strip()
    dataset_digest, dataset_bytes = dataset_identity(args.dataset)
    record = {
        "candidate": args.candidate,
        "configuration": configuration,
        "protocol": registry["protocol"],
        "dataset": str(args.dataset.resolve()),
        "dataset_sha256": dataset_digest,
        "dataset_bytes": dataset_bytes,
        "git_revision": revision,
        "elapsed_seconds": time.time() - started,
        **selection_score(output, registry["protocol"]),
    }
    (output / "experiment.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
