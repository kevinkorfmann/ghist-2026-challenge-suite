"""Versioned input guards for Task 1 calibration; not scientific certification."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "bottleneck_calibration_inputs_v1"
PARAMETERS = ("generations", "post_decline_fraction", "ancestral_effective_size")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_manifest(producer):
    paths = {
        "scripts/bottleneck_consumer_guard.py", producer,
        "src/ghist_bottleneck/provenance.py",
        "src/ghist_bottleneck/optimization.py",
        "src/ghist_bottleneck/sfs_baseline.py",
        "src/ghist_bottleneck/dadi_reference.py",
        "src/ghist_inference/validation.py",
    }
    if producer in {"scripts/evaluate_moments_chromosomes.py", "scripts/evaluate_bottleneck_ensemble.py"}:
        paths.add("scripts/evaluate_bottleneck_chromosomes.py")
    if producer == "scripts/evaluate_multiscale_sfs_recovery.py":
        paths.add("scripts/run_multiscale_sfs_recovery.py")
    return {name: sha256(ROOT / name) for name in sorted(paths)}


def capture_inputs(*paths):
    return {str(Path(path).resolve()): sha256(path) for path in paths}


def require_unchanged(hashes):
    for name, expected in hashes.items():
        if sha256(name) != expected:
            raise ValueError(f"Input or source changed during evaluation: {name}")


def require_new_outputs(*paths):
    if len({Path(path).resolve() for path in paths}) != len(paths):
        raise ValueError("Calibration output paths must be distinct")
    if any(Path(path).exists() for path in paths):
        raise FileExistsError("Refusing to overwrite existing calibration evidence")


def scalar_json(array, key):
    if key not in array or np.asarray(array[key]).shape != ():
        raise ValueError(f"Missing versioned calibration provenance: {key}")
    return json.loads(str(array[key]))


def physical(fit):
    values = np.asarray([fit[name] for name in PARAMETERS], dtype=float)
    if (not np.isfinite(values).all() or np.any(values <= 0)
            or not 0 < values[1] < 1):
        raise ValueError("Fit parameters must be finite and physically valid")
    return values


def latent(fit):
    generation, fraction, _ = physical(fit)
    return np.asarray([np.log(generation), np.log(fraction / (1 - fraction))])


def fit_data_contract(fit):
    values = np.asarray([fit.get("mutation_rate", np.nan), fit.get("sequence_length", np.nan)], dtype=float)
    if fit.get("folded") is not True or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("Calibration requires folded fits with verified positive mutation rate and exposure")
    if fit.get("callable_sequence_length", values[1]) != values[1]:
        raise ValueError("Complete-exposure calibration cannot be applied to a masked observation")
    sample = fit.get("sample_haplotypes")
    if type(sample) is not int or sample < 2:
        raise ValueError("Calibration requires a verified haplotype projection/sample size")
    if "sfs" in fit and fit["sfs"].get("projection") != sample:
        raise ValueError("Observed VCF projection differs from fitted sample size")
    return {"folded": True, "mutation_rate": float(values[0]), "sequence_length": float(values[1]),
            "sample_haplotypes": sample}


def require_fit(fit, *, engine=None, expected=None, expected_data=None):
    from ghist_bottleneck.optimization import validate_optimizer_contract

    contract = validate_optimizer_contract(
        fit, engine=engine, require_converged=True,
        require_agreement=True, require_neighborhood=True,
    )
    if expected is not None and contract != expected:
        raise ValueError("Calibration and observed optimizer contracts differ")
    data_contract = fit_data_contract(fit)
    if expected_data is not None and data_contract != expected_data:
        raise ValueError("Calibration and observed folding/mutation-rate/exposure/projection contracts differ")
    physical(fit)
    return contract


def bundle_metadata(fits, *, producer, inputs=None, simulation_design=None):
    if not fits:
        raise ValueError("Calibration needs all planned fits")
    contract = require_fit(fits[0])
    data_contract = fit_data_contract(fits[0])
    for fit in fits[1:]:
        require_fit(fit, expected=contract, expected_data=data_contract)
    return {
        "consumer_schema": np.asarray(SCHEMA),
        "optimizer_contract_json": np.asarray(json.dumps(contract, sort_keys=True)),
        "fit_records_json": np.asarray(json.dumps(fits, allow_nan=False)),
        "artifact_context_json": np.asarray(json.dumps({
            "producer": producer, "source_sha256": source_manifest(producer),
            "inputs": inputs or {}, "fit_data_contract": data_contract,
            "simulation_design": simulation_design,
        }, sort_keys=True)),
    }


def require_bundle(array, *, producer, engine=None, estimate_field=None, latent_field=None):
    if str(array.get("consumer_schema", "")) != SCHEMA:
        raise ValueError("Legacy or unknown calibration artifact; regenerate with current optimizer")
    context = scalar_json(array, "artifact_context_json")
    if (context.get("producer") != producer
            or context.get("source_sha256") != source_manifest(producer)):
        raise ValueError("Calibration producer source identity mismatch")
    contract = scalar_json(array, "optimizer_contract_json")
    fits = scalar_json(array, "fit_records_json")
    ids = unique_ids(array["cell_id"], "calibration cells")
    if not isinstance(fits, list) or len(fits) != len(ids):
        raise ValueError("Calibration fit records do not cover all cells")
    if not isinstance(context.get("fit_data_contract"), dict):
        raise ValueError("Calibration is missing its folding/rate/exposure contract")
    for index, fit in enumerate(fits):
        require_fit(fit, engine=engine, expected=contract, expected_data=context["fit_data_contract"])
        for field, expected in ((estimate_field, physical(fit)), (latent_field, latent(fit))):
            if field is not None:
                values = np.asarray(array[field], dtype=float)
                if (values.shape != (len(ids), len(expected))
                        or not np.allclose(values[index], expected, rtol=1e-12, atol=1e-12)):
                    raise ValueError(f"Calibration array differs from retained fit: {field}")
    return contract, fits, context


def unique_ids(values, context):
    ids = np.asarray(values)
    if (ids.ndim != 1 or ids.dtype.kind not in "iu" or not len(ids)
            or np.any(ids < 0) or len(np.unique(ids)) != len(ids)):
        raise ValueError(f"{context} must contain unique nonnegative integer cell IDs")
    return ids


def require_partition(split_record, all_cells, *, required_cells=None):
    all_cells = unique_ids(all_cells, "complete input")
    splits = split_record.get("split_cells", {})
    if set(splits) != {"calibration", "selection", "evaluation"}:
        raise ValueError("All three required data splits must be present")
    validated = {name: unique_ids(values, f"{name} split") for name, values in splits.items()}
    joined = np.concatenate(list(validated.values()))
    if len(np.unique(joined)) != len(joined) or set(joined) != set(all_cells):
        raise ValueError("Split cells must partition the complete input exactly")
    required = set(all_cells if required_cells is None else unique_ids(required_cells, "required cells"))
    if not required <= set(all_cells):
        raise ValueError("Required cells are missing from the complete input")
    selected = {name: np.asarray([cell for cell in values if cell in required], dtype=np.int64)
                for name, values in validated.items()}
    if any(not len(values) for values in selected.values()):
        raise ValueError("Every required split must contain at least one cell")
    return selected


def require_learned_binding(path, split_record):
    if split_record.get("prediction_npz_sha256") != sha256(path):
        raise ValueError("Learned predictions are not hash-bound to the split record")
    if split_record.get("evaluator_sha256") != sha256(ROOT / "scripts/evaluate_bottleneck_chromosomes.py"):
        raise ValueError("Learned prediction evaluator source identity mismatch")


def require_ensemble_inputs(mechanistic, learned, split_record, observed):
    contract, _, context = require_bundle(
        mechanistic, producer="scripts/evaluate_moments_chromosomes.py",
        engine="moments", latent_field="raw_latent",
    )
    require_fit(observed, engine="moments", expected=contract,
                expected_data=context["fit_data_contract"])
    learned_observed = split_record.get("observed_summary", {})
    physical({"generations": learned_observed.get("generations", np.nan),
              "post_decline_fraction": learned_observed.get("post_decline_fraction", np.nan),
              "ancestral_effective_size": 1.0})
    learned_ids = unique_ids(learned["cell_id"], "learned cells")
    scenarios = np.asarray(learned["scenario"])
    if scenarios.shape != learned_ids.shape:
        raise ValueError("Learned scenarios must cover every cell")
    required = learned_ids[scenarios == 0]
    splits = require_partition(split_record, learned_ids, required_cells=required)
    mechanical_ids = unique_ids(mechanistic["cell_id"], "mechanistic cells")
    if set(mechanical_ids) != set(required):
        raise ValueError("Mechanistic inputs must contain every nominal cell exactly once")
    learned_index = {int(cell): index for index, cell in enumerate(learned_ids)}
    mechanical_index = {int(cell): index for index, cell in enumerate(mechanical_ids)}
    for field in ("target_latent", "summary_calibrated_latent"):
        values = np.asarray(learned[field], dtype=float)
        if values.shape != (len(learned_ids), 2) or not np.isfinite(values).all():
            raise ValueError(f"Learned {field} must cover all cells with finite values")
    target = np.asarray(mechanistic["target_latent"], dtype=float)
    if target.shape != (len(mechanical_ids), 2) or not np.isfinite(target).all():
        raise ValueError("Mechanistic targets must cover all cells with finite values")
    for cell in mechanical_ids:
        if not np.array_equal(target[mechanical_index[int(cell)]],
                              learned["target_latent"][learned_index[int(cell)]]):
            raise ValueError("Ensemble inputs disagree on the target for a cell")
    return contract, context, splits
