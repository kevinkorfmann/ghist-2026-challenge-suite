"""Small atomic checkpoints with explicit configuration and content checks."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np


def atomic_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".")
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(record, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def spectrum_digest(spectrum):
    values = np.asarray(spectrum, dtype="<f8", order="C")
    return hashlib.sha256(str(values.shape).encode() + values.tobytes()).hexdigest()


def save_spectrum_checkpoint(path, spectrum, metadata):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(metadata, spectrum_sha256=spectrum_digest(spectrum))
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".")
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(
                handle, spectrum=spectrum, metadata=json.dumps(metadata)
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return metadata


def load_spectrum_checkpoint(path, *, expected):
    with np.load(path, allow_pickle=False) as record:
        spectrum = np.asarray(record["spectrum"], dtype=float)
        metadata = json.loads(str(record["metadata"]))
    for key, value in expected.items():
        # JSON normalizes tuples to lists and numpy-free scalar configuration.
        if metadata.get(key) != json.loads(json.dumps(value)):
            raise ValueError(f"Spectrum checkpoint configuration mismatch: {key}")
    if metadata.get("spectrum_sha256") != spectrum_digest(spectrum):
        raise ValueError("Spectrum checkpoint digest mismatch")
    return spectrum, metadata


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def simulation_identity(arguments, *, stdpopsim_version, scaling_scheme):
    excluded = {
        "out",
        "trees",
        "coding_bed",
        "slim_path",
        "fit_restarts",
        "skip_fit",
        "resume",
        "existing_trees_sha256",
    }
    return {
        "parameters": {
            key: value for key, value in vars(arguments).items() if key not in excluded
        },
        "coding_bed_sha256": file_digest(arguments.coding_bed),
        "stdpopsim_version": stdpopsim_version,
        "simulation_scaling_scheme": scaling_scheme,
    }


def simulation_tree_checkpoint(arguments, identity, simulate):
    """Resume a verified tree, or record a new simulation before inference.

    Importing a pre-checkpoint tree requires its externally recorded SHA-256.
    The caller supplies the original simulation arguments; the import marker
    makes that provenance distinction explicit instead of inventing metadata.
    """
    import tskit

    manifest = arguments.out.with_suffix(".simulation.json")
    trees = arguments.trees
    if arguments.resume:
        digest = file_digest(trees)
        if manifest.exists():
            record = json.loads(manifest.read_text())
            if record["simulation"] != identity or record["tree_sha256"] != digest:
                raise ValueError(
                    "Simulation checkpoint configuration or tree digest mismatch"
                )
        else:
            if arguments.existing_trees_sha256 != digest:
                raise ValueError(
                    "Importing an existing tree requires its recorded SHA-256"
                )
            atomic_json(
                manifest,
                {
                    "simulation": identity,
                    "tree_sha256": digest,
                    "imported_existing_tree": True,
                    "provenance": "Caller supplied original simulation arguments and externally recorded tree SHA-256",
                },
            )
        return tskit.load(trees)
    if trees.exists() or manifest.exists():
        raise FileExistsError(
            "Simulation output already exists; use --resume after checking its configuration"
        )
    result = simulate()
    trees.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=trees.parent, prefix=trees.name + ".")
    os.close(descriptor)
    try:
        result.dump(temporary)
        os.replace(temporary, trees)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    atomic_json(
        manifest,
        {
            "simulation": identity,
            "tree_sha256": file_digest(trees),
            "imported_existing_tree": False,
        },
    )
    return result


def fit_checkpoint(path, *, spectrum, fit, settings, resume):
    specification = {
        "protocol": "bounded_refinement_diagnostics_v1",
        "spectrum_sha256": spectrum_digest(spectrum),
        "settings": json.loads(json.dumps(settings)),
    }
    path = Path(path)
    if path.exists():
        record = json.loads(path.read_text())
        if not resume or record["specification"] != specification:
            raise ValueError(
                "Fit checkpoint exists with incompatible settings or without --resume"
            )
        if record["status"] == "completed":
            return record["fit"]
    try:
        result = fit(spectrum, **settings)
    except Exception as error:
        atomic_json(
            path,
            {
                "specification": specification,
                "status": "failed",
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise
    atomic_json(
        path, {"specification": specification, "status": "completed", "fit": result}
    )
    return result
