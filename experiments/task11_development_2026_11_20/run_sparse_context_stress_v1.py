"""Frozen-model sparse-family stress; synthetic inputs only, no refitting."""
import argparse
import hashlib
import importlib.util
import json
import os
from collections import Counter
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def verify_hashes(bound):
    for path, expected in bound.items():
        if sha(path) != expected:
            raise ValueError(f"Bound input changed: {path}")


def initialize(root):
    candidate_path = root / "sparse_context_candidate_v1.json"
    candidate = load(candidate_path)
    if candidate["schema"] != "task11-frozen-expanded-context-stress-v1" or not candidate["frozen_before_stress"]:
        raise ValueError("Unknown or unfrozen candidate")
    if candidate["runner_sha256"] != sha(__file__):
        raise ValueError("Stress runner changed")
    verify_hashes(candidate["bound_inputs"])
    guard = load(root / "expanded_development_audit_v1/audit.json")
    if not guard["sparse_stress_launch_allowed"]:
        raise ValueError("Expanded development gates did not pass")
    protocol = load(root / "sparse_context_generation_protocol_v1.json")
    design = load(root / "sparse_context_stress_design_v1.json")
    if protocol["cohorts"] != design["cohorts"] or protocol["gates"] != design["gates"]:
        raise ValueError("Stress generation differs from prospective design")
    report = load(root / "expanded_context_fit_v2/report.json")
    for role in ("base", "context"):
        if candidate[f"{role}_temperature"] != report[f"{role}_temperature"]:
            raise ValueError("Frozen calibrated temperature changed")
        model_path = Path(candidate[f"{role}_model"]).resolve()
        if model_path != root / "expanded_context_fit_v2" / f"{role}_model.joblib":
            raise ValueError("Unexpected model path")
        if candidate["bound_inputs"].get(str(model_path)) != report["outputs_sha256"][model_path.name]:
            raise ValueError("Model does not match completed expanded fit")
    source = root / "source_v5"
    verify_hashes({str(source / path): value for path, value in protocol["source_sha256"].items()})
    helper = root / "context_diagnostic_source_v1/diagnose_family_context.py"
    spec = importlib.util.spec_from_file_location("frozen_context", helper)
    ctx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ctx)
    dev = ctx.dev
    if dev.ROOT.resolve() != source.resolve():
        raise ValueError("Unexpected source module")
    dev.validate_protocol(protocol)
    dev.require_runtime(protocol)
    dev.require_hashes(protocol["inputs"])
    if "SLURM_JOB_ID" not in os.environ:
        raise RuntimeError("Stress compute requires Slurm")
    return candidate, protocol, ctx


def infer(directory, candidate, ctx, map_path):
    """Only observed feature pairs, detector segments and frozen models enter inference."""
    dev, np = ctx.dev, ctx.np
    pairs = dev.load(directory / "feature_pairs.json")
    features = np.load(directory / "features.npz")["features"]
    if features.shape != (len(pairs), 149):
        raise ValueError("Expected all 149 frozen v5 pair features")
    base_model = dev.joblib.load(candidate["base_model"])
    context_model = dev.joblib.load(candidate["context_model"])
    base = dev.predict(pairs, features, base_model, candidate["base_temperature"], False)
    samples = {sample for row in pairs for sample in dev.key(row)}
    records, _ = dev.read_haplotype_segments(directory / "segments.ibd.gz", samples)
    masks = {cutoff: dev.read_hbd_masks(directory / "segments.hbd.gz", cutoff)
             for cutoff in ctx.PARAMETERS["hbd_cutoffs_cm"]}
    extra, audit = ctx.context_features(pairs, ctx.probabilities(base), list(dev.LABELS).index("siblings"),
                                        records, masks, dev.read_map(map_path))
    context = dev.predict(pairs, np.concatenate([features, extra], axis=1), context_model,
                          candidate["context_temperature"], False)
    return base, context, extra, audit


def run(root, index, candidate, protocol, ctx):
    dev = ctx.dev
    config = next(c for c in protocol["cohorts"] if c["index"] == index)
    if config["split"] != "sparse_stress":
        raise ValueError("Only frozen synthetic sparse stress is allowed")
    project = root.parents[1]
    map_path = project / "external/dog-recombination/campbell2016_average.plink.map"
    directory = root / "results_sparse_context_v1" / f"cohort_{index:02d}"
    result = dev.run_cohort(argparse.Namespace(protocol=root / "sparse_context_generation_protocol_v1.json",
        split="sparse_stress", index=index, map=map_path, hapibd=project / "external/hap-ibd/hap-ibd.jar",
        baseline=None, candidate=None, out=directory))
    if not result["completed"]:
        raise RuntimeError(f"Sparse generation failed: {directory}")
    output = root / "predictions_sparse_context_v1" / f"cohort_{index:02d}"
    output.mkdir(parents=True, exist_ok=False)
    base, context, extra, grouping = infer(directory, candidate, ctx, map_path)
    dev.dump(output / "base_predictions.json", base)
    dev.dump(output / "context_predictions.json", context)
    dev.dump(output / "grouping_audit.json", grouping)
    ctx.np.savez_compressed(output / "context_features.npz", features=extra)
    truth = dev.load(directory / "truth.json")
    verify_hashes(candidate["bound_inputs"])
    dev.require_hashes({str(directory / path): value for path, value in result["output_sha256"].items()})
    audit = {"cohort": config, "candidate_sha256": sha(root / "sparse_context_candidate_v1.json"),
        "generation_result_sha256": sha(directory / "result.json"),
        "base_score": dev.score(truth, base), "context_score": dev.score(truth, context),
        "output_sha256": {path.name: sha(path) for path in sorted(output.iterdir())},
        "truth_used_for_inference": False, "slurm_job_id": os.environ["SLURM_JOB_ID"], "submission_ready": False}
    dev.dump(output / "audit.json", audit)
    print(json.dumps(audit, indent=2))


def summarize(root, candidate, protocol, ctx):
    dev, np = ctx.dev, ctx.np
    gates = dev.validate_protocol(protocol)
    out = root / "sparse_context_summary_v1"
    out.mkdir(exist_ok=False)
    cohorts, partitions, bounds = [], [], {}
    for config in protocol["cohorts"]:
        index = config["index"]
        directory = root / "results_sparse_context_v1" / f"cohort_{index:02d}"
        prediction = root / "predictions_sparse_context_v1" / f"cohort_{index:02d}"
        data = dev.dataset(directory / "result.json")
        if data[0]["cohort"] != config:
            raise ValueError("Unexpected stress cohort")
        master = dev.load(directory / "protocol.json")
        if master["master_sha256"] != sha(root / "sparse_context_generation_protocol_v1.json"):
            raise ValueError("Generation protocol mismatch")
        audit = dev.load(prediction / "audit.json")
        if audit["candidate_sha256"] != sha(root / "sparse_context_candidate_v1.json"):
            raise ValueError("Candidate mismatch")
        if audit["generation_result_sha256"] != sha(directory / "result.json"):
            raise ValueError("Generation result changed")
        verify_hashes({str(prediction / path): value for path, value in audit["output_sha256"].items()})
        bounds[str(prediction / "audit.json")] = sha(prediction / "audit.json")
        bounds[str(directory / "result.json")] = sha(directory / "result.json")
        rows = dev.load(prediction / "context_predictions.json")
        if dev.score(data[1], rows) != audit["context_score"]:
            raise ValueError("Stored stress score mismatch")
        cohorts.append({"cohort_id": f"sparse_stress_{index:02d}", "truth": data[1]["pairs"], "predictions": rows})
        partitions.append(dev.partition_metrics(data[1]["families"], dev.load(directory / "clusters.json")["families"]))
    truth, rows, per_cohort = {"pairs": []}, [], []
    for index, cohort in enumerate(cohorts):
        per_cohort.append(dev.score({"pairs": cohort["truth"]}, cohort["predictions"]))
        for source, target in ((cohort["truth"], truth["pairs"]), (cohort["predictions"], rows)):
            target.extend({**row, "indID1": f"c{index}_{row['indID1']}", "indID2": f"c{index}_{row['indID2']}"} for row in source)
    overall = dev.score(truth, rows)
    rng = np.random.default_rng(20301120)
    bootstrap = []
    for _ in range(2000):
        correct, counts = Counter(), Counter()
        for index in rng.integers(0, len(cohorts), len(cohorts)):
            correct.update(per_cohort[index]["correct_counts"])
            counts.update(per_cohort[index]["counts"])
        bootstrap.append(np.mean([correct[label] / counts[label] for label in dev.LABELS]))
    overall["cohort_bootstrap_balanced_interval95"] = np.quantile(bootstrap, [.025, .975]).tolist()
    exact = sum(p["exactly_recovered_families"] for p in partitions) / (8 * len(cohorts))
    missing = sum(len(p["missing_samples"]) for p in partitions)
    checks = {"accuracy": overall["accuracy"] >= gates["accuracy_min"],
        "balanced_accuracy": overall["balanced_accuracy"] >= gates["balanced_accuracy_min"],
        "each_label_recall": min(overall["recall"].values()) >= gates["each_label_recall_min"],
        "bootstrap_lower": overall["cohort_bootstrap_balanced_interval95"][0] >= gates["bootstrap_balanced_lower95_min"],
        "family_recovery": exact >= gates["exact_family_recovery_min"],
        "missing_samples": missing <= gates["missing_samples_max"],
        "missing_true_pairs": overall["missing_true_pairs"] <= gates["missing_true_pairs_max"],
        "extra_pairs": overall["extra_cross_family_pairs"] <= gates["extra_cross_family_pairs_max"],
        "calibration": overall["ece_10_equal_width"] <= gates["expected_calibration_error_max"]}
    dev.dump(out / "independent_scoring_input.json", {"labels": list(dev.LABELS), "cohorts": cohorts})
    verify_hashes(candidate["bound_inputs"])
    verify_hashes(bounds)
    summary = {"scope": "Required sparse-context stress, separate from reserved confirmation", "score": overall,
        "checks": checks, "all_original_gates_passed": all(checks.values()), "exact_family_recovery": exact,
        "candidate_sha256": sha(root / "sparse_context_candidate_v1.json"), "bound_results": bounds,
        "scoring_input_sha256": sha(out / "independent_scoring_input.json"), "cohorts": protocol["cohorts"],
        "submission_ready": False, "fresh_confirmation_required": True}
    dev.dump(out / "summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("command", choices=("run", "summarize"))
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    root = args.root.resolve()
    candidate, protocol, ctx = initialize(root)
    if args.command == "run":
        run(root, args.index, candidate, protocol, ctx)
    else:
        summarize(root, candidate, protocol, ctx)


if __name__ == "__main__":
    main()
