"""Development-only context-wrapper parity on one retained calibration cohort."""
import importlib.util
import json
import os
import sys
from pathlib import Path


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


root = Path(sys.argv[1]).resolve()
wrapper = module("stress_wrapper", root / "run_sparse_context_stress_v1.py")
helper = root / "context_diagnostic_source_v1/diagnose_family_context.py"
if wrapper.sha(helper) != "4b347c189e0a5799c99e158715c236934485e618db44cb231e84a27535fd396b":
    raise ValueError("Context helper changed")
ctx = module("frozen_context", helper)
dev = ctx.dev
if "SLURM_JOB_ID" not in os.environ:
    raise RuntimeError("Inference requires Slurm")
protocol = dev.load(root / "protocol_v5.json")
dev.require_runtime(protocol)
dev.require_hashes({str(dev.ROOT / path): value for path, value in protocol["source_sha256"].items()})
report = dev.load(root / "expanded_context_fit_v2/report.json")
bound = {str(root / "expanded_context_fit_v2" / path): value for path, value in report["outputs_sha256"].items()}
dev.require_hashes(bound)
directory = root / "results_v5/calibration_00"
result = dev.load(directory / "result.json")
dev.require_hashes({str(directory / path): value for path, value in result["output_sha256"].items()})
candidate = {f"{role}_model": str(root / f"expanded_context_fit_v2/{role}_model.joblib") for role in ("base", "context")}
candidate.update({f"{role}_temperature": report[f"{role}_temperature"] for role in ("base", "context")})
base, context, features, grouping = wrapper.infer(directory, candidate, ctx,
    root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map")
checks = {}
for role, predictions in (("base", base), ("context", context)):
    retained = dev.load(root / f"expanded_context_fit_v2/scoring_calibration_{role}.json")["cohorts"][0]["predictions"]
    expected = {dev.key(row): row for row in retained}
    actual = {dev.key(row): row for row in predictions}
    if actual.keys() != expected.keys():
        raise ValueError("Wrapper pair coverage differs")
    changes = sum(actual[pair]["relation"] != expected[pair]["relation"] for pair in expected)
    delta = max(abs(actual[pair]["probabilities"][label] - expected[pair]["probabilities"][label])
                for pair in expected for label in dev.LABELS)
    checks[role] = {"pairs": len(expected), "label_changes": changes, "maximum_probability_delta": delta}
    if changes or delta > 1e-12:
        raise ValueError("Wrapper changed frozen predictions")
out = root / "expanded_context_smoke_v1"
out.mkdir(exist_ok=False)
dev.dump(out / "base_predictions.json", base)
dev.dump(out / "context_predictions.json", context)
dev.dump(out / "grouping_audit.json", grouping)
ctx.np.savez_compressed(out / "context_features.npz", features=features)
dev.require_hashes(bound)
audit = {"scope": "One reused calibration cohort; model inference parity, not new scientific evidence",
    "checks": checks, "all_passed": True, "detector_rerun": False, "truth_used_for_inference": False,
    "report_sha256": dev.digest(root / "expanded_context_fit_v2/report.json"),
    "wrapper_sha256": dev.digest(root / "run_sparse_context_stress_v1.py"), "smoke_script_sha256": dev.digest(__file__),
    "output_sha256": {path.name: dev.digest(path) for path in sorted(out.iterdir())}, "submission_ready": False}
dev.dump(out / "audit.json", audit)
print(json.dumps(audit, indent=2))
