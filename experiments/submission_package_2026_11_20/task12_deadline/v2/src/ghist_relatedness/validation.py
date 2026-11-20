"""Provenance guards for independent relatedness simulation splits."""

from __future__ import annotations

import hashlib
from pathlib import Path
from collections.abc import Mapping, Sequence


def validate_simulation_splits(
    splits: Mapping[str, Sequence[Path]],
) -> dict[str, list[dict[str, str]]]:
    """Reject repeated files, including renamed copies, before fitting.

    Distinct hashes cannot establish independent seeds or pedigrees by themselves;
    callers must also preserve simulation provenance and grouped split design.
    """
    seen: dict[str, tuple[str, Path]] = {}
    manifest: dict[str, list[dict[str, str]]] = {}
    for role, paths in splits.items():
        if not paths:
            raise ValueError(f"Simulation split {role} is empty")
        manifest[role] = []
        for raw_path in paths:
            path = Path(raw_path)
            if not path.stat().st_size:
                raise ValueError(f"Simulation input is empty: {path}")
            with path.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            if digest in seen:
                previous_role, previous_path = seen[digest]
                raise ValueError(
                    f"Duplicate simulation content in {previous_role} ({previous_path}) and {role} ({path}); independent splits are required"
                )
            seen[digest] = (role, path)
            manifest[role].append({"path": str(path), "sha256": digest})
    return manifest
