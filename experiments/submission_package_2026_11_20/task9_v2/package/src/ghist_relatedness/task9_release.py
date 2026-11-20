"""Explicit, checksum-pinned Task 9 release and input contracts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import sklearn

from ghist_relatedness.labels import IBD_LABEL_MEIOSES


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def input_pairs(rows):
    """Validate raw segments and derive coverage before feature extraction."""
    pairs, seen = set(), set()
    for row in rows:
        pair = tuple(sorted((row["indID1"], row["indID2"])))
        if any(not name.strip() or name != name.strip() for name in pair) or pair[0] == pair[1]:
            raise ValueError("Segments require distinct, nonempty, unpadded sample IDs")
        if row["type"] not in {"IBD1", "IBD2"}:
            raise ValueError("Task 9 permits IBD1 and IBD2 segments only")
        numbers = [float(row[name]) for name in ("start_bp", "end_bp", "start_cM", "end_cM", "length")]
        start_bp, end_bp, start_cm, end_cm, length = numbers
        if (
            not np.all(np.isfinite(numbers))
            or min(numbers) < 0
            or end_bp < start_bp
            or end_cm < start_cm
            or not np.isclose(end_cm - start_cm, length, rtol=1e-5, atol=1e-4)
        ):
            raise ValueError("Invalid or inconsistent segment coordinates/length")
        key = (*pair, row["chr"], row["type"], *numbers)
        if key in seen:
            raise ValueError("Duplicate IBD segment, including reversed pair order")
        seen.add(key)
        pairs.add(pair)
    if not pairs:
        raise ValueError("No IBD pairs in input")
    return pairs


def load_release(path: Path, artifact_root: Path, sex_map: Path, genome_cm=None):
    manifest = json.loads(path.read_text())
    if (
        manifest.get("schema") != "ghist.task9.release.v1"
        or manifest.get("status") != "VALIDATED_RELEASE_CANDIDATE"
        or set(manifest.get("accepted_labels", [])) != set(IBD_LABEL_MEIOSES)
        or manifest.get("pair_domain") != "unordered_pairs_present_in_segment_csv"
        or not manifest.get("gates")
        or not all(value is True for value in manifest["gates"].values())
    ):
        raise ValueError("Invalid or incomplete Task 9 release manifest")
    if manifest["environment"]["scikit_learn"] != sklearn.__version__:
        raise ValueError("Task 9 release requires its recorded scikit-learn version")
    model_path = artifact_root / manifest["model"]["path"]
    if sha256(model_path) != manifest["model"]["sha256"]:
        raise ValueError("Task 9 release model checksum mismatch")
    if sha256(sex_map) != manifest["sex_map_sha256"]:
        raise ValueError("Task 9 release genetic-map checksum mismatch")
    for name in ("ibd.py", "labels.py"):
        if sha256(Path(__file__).with_name(name)) != manifest["feature_contract_sha256"][name]:
            raise ValueError(f"Task 9 release feature contract changed: {name}")
    expected_cm = float(manifest["genome_cm"])
    if genome_cm is not None and genome_cm != expected_cm:
        raise ValueError("Task 9 release genome normalization mismatch")
    return manifest, model_path, expected_cm


def validate_payload(payload, manifest=None):
    model, features = payload["model"], payload["features"]
    if len(features) != len(set(features)) or len(features) != model.n_features_in_:
        raise ValueError("Invalid model feature contract")
    if manifest is not None:
        if set(model.classes_) != set(IBD_LABEL_MEIOSES):
            raise ValueError("Task 9 release must contain exactly the accepted labels")
        if list(features) != manifest["features"]:
            raise ValueError("Task 9 release feature order mismatch")
    return model, features
