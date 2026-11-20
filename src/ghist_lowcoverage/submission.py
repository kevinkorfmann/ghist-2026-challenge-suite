"""Exact official-schema formatting for review drafts, without uploading."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import yaml


TEMPLATE_SHA256 = {
    "bottleneck": "d2d4ad595d4de0466fc4291d3f6550054da98676e42c662b6b57aaeabf5a9936",
    "split_isolation": "46bddd26dd3ef9fac2d707c549214b5b8f9d2c427b9d0a38cc7ecc9a4ed2d0cb",
}
PARAMETERS = {
    "bottleneck": {"generations", "post_decline_fraction"},
    "split_isolation": {"pop1_population_size", "pop2_population_size", "generations_since_split"},
}


def format_submission_draft(task: str, targets: dict, template: Path) -> str:
    """Fill a verified template; formatting is not scientific authorization."""
    if task not in PARAMETERS or set(targets) != PARAMETERS[task]:
        raise ValueError("unknown task or parameter schema differs from official template")
    data = template.read_bytes()
    if hashlib.sha256(data).hexdigest() != TEMPLATE_SHA256[task]:
        raise ValueError("official template checksum mismatch")
    record = yaml.safe_load(data)
    if set(record["parameters"]) != PARAMETERS[task]:
        raise ValueError("template parameter schema mismatch")
    for name, value in targets.items():
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError("parameters must be finite positive numbers")
        if name == "post_decline_fraction" and value >= 1:
            raise ValueError("post-decline fraction must be less than one")
        record["parameters"][name] = float(value)
    return yaml.safe_dump(record, sort_keys=False)
