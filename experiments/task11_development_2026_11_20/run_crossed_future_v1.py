"""Reserved synthetic stress/confirmation with the fixed crossed context model.

The candidate is frozen only after all development launch gates pass. Nothing
here refits a model, selects a threshold, accesses final data, or submits.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import traceback
from collections import Counter
from pathlib import Path

PLAN_SHA = "e90d2235de805d17a6c1d0176cab2132c2626b5ef04d3c97b4a4e12c67850016"
ADAPTER_SHA = "994b268b117e047fcfbc42eef54b5dcbbe56da2602e043746f5200346bbcbbe4"
INFER_SHA = "49bfb17e89352669f74b2ad72e7256251482674c147f9f19f8729155b8824a56"
CHECKS = {"accuracy", "balanced_accuracy", "each_label_recall", "bootstrap_lower",
          "family_recovery", "missing_samples", "missing_true_pairs", "extra_pairs", "calibration"}


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify(bound):
    if not isinstance(bound, dict) or not bound:
        raise ValueError("Missing hash bindings")
    for path, expected in bound.items():
        if sha(path) != expected:
            raise ValueError(f"Bound input changed: {path}")


def import_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_development(audit):
    """Keep the exact prospective launch guard, including its time budget."""
    dev = audit["summaries"]["development_context"]
    cal = audit["summaries"]["calibration_context"]["checks"]
    if (audit.get("input_design_sha256") != PLAN_SHA or
            audit.get("fresh_stress_launch_allowed") is not True or
            audit.get("scientific_development_gates_passed") is not True or
            audit.get("within_decision_budget") is not True or
            set(dev["checks"]) != CHECKS or any(v is not True for v in dev["checks"].values()) or
            dev["all_nine_gates_passed"] is not True or
            cal.get("each_label_recall") is not True or cal.get("calibration") is not True):
        raise ValueError("The fixed crossed development/calibration launch gates failed")


def require_stress(summary, candidate_sha):
    if (summary.get("split") != "stress" or summary.get("candidate_sha256") != candidate_sha or
            summary.get("plan_sha256") != PLAN_SHA or summary.get("all_original_gates_passed") is not True or
            set(summary.get("checks", {})) != CHECKS or
            any(v is not True for v in summary["checks"].values()) or
            len(summary.get("cohorts", [])) != 8):
        raise ValueError("Passing candidate-bound stress evidence is required before confirmation")


def initialize(root, candidate_path=None, candidate_sha=None):
    """Check the source closure before importing any scientific candidate code."""
    plan_path, protocol_path = root / "crossed_development_plan_v1.json", root / "protocol_v5.json"
    verify({str(plan_path): PLAN_SHA, str(root / "adapter_dependency_v3.py"): ADAPTER_SHA,
            str(root / "run_sparse_context_stress_v1.py"): INFER_SHA})
    plan, protocol = load(plan_path), load(protocol_path)
    if plan["schema"] != "task11-crossed-development-v1" or plan["gates"] != protocol["gates"]:
        raise ValueError("Changed plan/schema/gates")
    verify(plan["bound_inputs"])
    adapter = import_path("crossed_adapter_dependency", root / "adapter_dependency_v3.py")
    source = root / "source_v5"
    adapter.verify_source(source, protocol["source_sha256"])
    verify(protocol["inputs"])
    candidate = None
    if candidate_path is not None:
        verify({str(candidate_path): candidate_sha})
        candidate = load(candidate_path)
        if (candidate["schema"] != "task11-fixed-crossed-context-candidate-v1" or
                candidate["plan_sha256"] != PLAN_SHA or candidate["gates"] != plan["gates"] or
                candidate["frozen_before_stress"] is not True or
                candidate["bound_inputs"].get(str(Path(__file__).resolve())) != sha(__file__)):
            raise ValueError("Unexpected candidate identity or runner")
        verify(candidate["bound_inputs"])
        audit_path = root / "crossed_development_audit_v1/audit.json"
        require_development(load(audit_path))
        report = load(root / "crossed_context_fit_v1/report.json")
        for role in ("base", "context"):
            model = root / "crossed_context_fit_v1" / f"{role}_model.joblib"
            if (Path(candidate[f"{role}_model"]) != model or
                    candidate[f"{role}_temperature"] != report[f"{role}_temperature"] or
                    candidate["bound_inputs"].get(str(model)) != report["outputs_sha256"][model.name]):
                raise ValueError("Frozen model or calibrated temperature differs from fit")
    if any(name == "scripts" or name.startswith(("scripts.", "ghist_")) for name in sys.modules):
        raise ValueError("Scientific modules were imported before verification")
    sys.dont_write_bytecode = True
    paths = [str(source), str(source / "src"), str(source / "scripts")]
    sys.path[:0] = paths
    os.environ["PYTHONPATH"] = os.pathsep.join(paths)
    ctx = import_path("crossed_frozen_context", root / "context_diagnostic_source_v1/diagnose_family_context.py")
    # The legacy helper accepts root as argv[1]; this adapter also passes it
    # explicitly so imported use does not depend on the surrounding CLI.
    import importlib
    ctx.dev = importlib.import_module("scripts.develop_task11_pedigree")
    if ctx.dev.ROOT.resolve() != source.resolve() or ctx.PARAMETERS != plan["context_parameters"]:
        raise ValueError("Unexpected scientific source or context feature specification")
    ctx.dev.validate_protocol(protocol)
    ctx.dev.require_runtime(protocol)
    infer = import_path("crossed_frozen_inference", root / "run_sparse_context_stress_v1.py").infer
    return plan, protocol, candidate, ctx, adapter, infer


def freeze(root):
    plan, protocol, _, ctx, _, _ = initialize(root)
    audit_path = root / "crossed_development_audit_v1/audit.json"
    audit = load(audit_path)
    require_development(audit)
    fit = root / "crossed_context_fit_v1"
    report_path, report = fit / "report.json", load(fit / "report.json")
    verify(report["inputs_sha256"])
    verify({str(fit / name): value for name, value in report["outputs_sha256"].items()})
    if audit["fit_report_sha256"] != sha(report_path):
        raise ValueError("Development audit/fit mismatch")
    bounds = {**plan["bound_inputs"], str(root / "crossed_development_plan_v1.json"): PLAN_SHA,
              str(audit_path): sha(audit_path), str(report_path): sha(report_path),
              str(root / "adapter_dependency_v3.py"): ADAPTER_SHA,
              str(root / "run_sparse_context_stress_v1.py"): INFER_SHA,
              str(Path(__file__).resolve()): sha(__file__),
              str(root / "prepare_crossed_context_v1.py"): sha(root / "prepare_crossed_context_v1.py")}
    candidate = {"schema": "task11-fixed-crossed-context-candidate-v1", "plan_sha256": PLAN_SHA,
        "protocol_sha256": sha(root / "protocol_v5.json"), "gates": plan["gates"],
        "source_sha256": protocol["source_sha256"], "frozen_before_stress": True,
        "frozen_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "context_features": 185, "base_features": 149, "bound_inputs": bounds,
        "submission_ready": False, "no_final_challenge_access": True}
    for role in ("base", "context"):
        model = fit / f"{role}_model.joblib"
        candidate[f"{role}_model"] = str(model)
        candidate[f"{role}_temperature"] = report[f"{role}_temperature"]
        bounds[str(model)] = sha(model)
    ctx.dev.require_hashes(bounds)
    destination = root / "crossed_context_candidate_v1.json"
    with destination.open("x") as stream:
        stream.write(json.dumps(candidate, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"candidate": str(destination), "sha256": sha(destination)}))


def check_launch(root, split, plan, candidate_sha, stress_sha):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Synthetic compute requires Slurm")
    if datetime.datetime.now(datetime.UTC) > datetime.datetime.fromisoformat(plan["budget"]["latest_stress_confirmation_start_utc"]):
        raise RuntimeError("Hard future-cohort launch deadline passed")
    if split == "confirmation":
        stress = root / "crossed_stress_summary_v1/summary.json"
        verify({str(stress): stress_sha})
        value = load(stress)
        require_stress(value, candidate_sha)
        if value["cohorts"] != plan["reserved_future_cohorts"]["stress"]:
            raise ValueError("Stress evidence uses a different registered design")
        verify(value["bound_results"])


def run(root, split, index, candidate_path, candidate_sha, stress_sha):
    plan, protocol, candidate, ctx, _, infer = initialize(root, candidate_path, candidate_sha)
    check_launch(root, split, plan, candidate_sha, stress_sha)
    config = plan["reserved_future_cohorts"][split][index]
    if config["index"] != index or config["split"] != split:
        raise ValueError("Unexpected reserved cohort")
    generator = import_path("crossed_reserved_generator", root / "crossed_source_v1/generate_crossed_genotypes.py")
    out = root / f"results_crossed_{split}_v1" / f"cohort_{index:02d}"
    out.mkdir(parents=True, exist_ok=False)
    dev = ctx.dev
    dump(out / "protocol.json", {"master_sha256": PLAN_SHA, "candidate_sha256": candidate_sha,
        "cohort": config, "job_id": os.environ["SLURM_JOB_ID"], "stress_summary_sha256": stress_sha,
        "source_sha256": candidate["bound_inputs"], "final_challenge_access": False})
    result = {"cohort": config, "candidate_sha256": candidate_sha,
        "protocol_sha256": sha(out / "protocol.json"), "submission_ready": False}
    try:
        project = root.parents[1]
        map_path = project / "external/dog-recombination/campbell2016_average.plink.map"
        truth = generator.generate_genotypes(map_path, out, config)
        truth_hash = sha(out / "truth.json")
        vcf_hash = sha(out / "synthetic.vcf.gz")
        inference = dev.detect_cluster_extract(out, map_path, project / "external/hap-ibd/hap-ibd.jar",
                                              truth["samples"], protocol["clustering_policy"])
        base, rows, extra, grouping = infer(out, candidate, ctx, map_path)
        dump(out / "base_predictions.json", base)
        dump(out / "context_predictions.json", rows)
        dump(out / "grouping_audit.json", grouping)
        ctx.np.savez_compressed(out / "context_features.npz", features=extra)
        result.update(inference=inference, partition=dev.partition_metrics(truth["families"], inference["clusters"]["families"]),
                      base_score=dev.score(truth, base), context_score=dev.score(truth, rows), truth_used_for_inference=False)
        verify({**candidate["bound_inputs"], str(out / "truth.json"): truth_hash, str(out / "synthetic.vcf.gz"): vcf_hash})
        result["completed"] = True
    except Exception as error:  # noqa: BLE001 - preserve every failed registered outcome
        result.update(completed=False, failure={"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()})
    result["output_sha256"] = {str(path.relative_to(out)): sha(path) for path in sorted(out.rglob("*")) if path.is_file()}
    dump(out / "result.json", result)
    if not result["completed"]:
        raise RuntimeError(f"Registered cohort failed; retained at {out}")
    print(json.dumps({"split": split, "index": index, "completed": True}))


def gate_checks(overall, exact, missing, gates):
    return {"accuracy": overall["accuracy"] >= gates["accuracy_min"],
        "balanced_accuracy": overall["balanced_accuracy"] >= gates["balanced_accuracy_min"],
        "each_label_recall": min(overall["recall"].values()) >= gates["each_label_recall_min"],
        "bootstrap_lower": overall["cohort_bootstrap_balanced_interval95"][0] >= gates["bootstrap_balanced_lower95_min"],
        "family_recovery": exact >= gates["exact_family_recovery_min"],
        "missing_samples": missing <= gates["missing_samples_max"],
        "missing_true_pairs": overall["missing_true_pairs"] <= gates["missing_true_pairs_max"],
        "extra_pairs": overall["extra_cross_family_pairs"] <= gates["extra_cross_family_pairs_max"],
        "calibration": overall["ece_10_equal_width"] <= gates["expected_calibration_error_max"]}


def summarize(root, split, candidate_path, candidate_sha):
    plan, _, candidate, ctx, _, _ = initialize(root, candidate_path, candidate_sha)
    dev, np = ctx.dev, ctx.np
    configs = plan["reserved_future_cohorts"][split]
    out = root / f"crossed_{split}_summary_v1"
    out.mkdir(exist_ok=False)
    cohorts, partitions, bounds, base_cohorts = [], [], {}, []
    for config in configs:
        directory = root / f"results_crossed_{split}_v1" / f"cohort_{config['index']:02d}"
        result = load(directory / "result.json")
        if not result["completed"] or result["cohort"] != config or result["candidate_sha256"] != candidate_sha:
            raise ValueError("Failed or mismatched reserved cohort")
        verify({str(directory / path): value for path, value in result["output_sha256"].items()})
        master = load(directory / "protocol.json")
        if master["master_sha256"] != PLAN_SHA or master["candidate_sha256"] != candidate_sha or master["cohort"] != config:
            raise ValueError("Reserved cohort protocol changed")
        bounds[str(directory / "result.json")] = sha(directory / "result.json")
        bounds.update({str(directory / name): result["output_sha256"][name] for name in
                       ("truth.json", "clusters.json", "protocol.json", "base_predictions.json", "context_predictions.json")})
        truth, clusters = load(directory / "truth.json"), load(directory / "clusters.json")
        partition = dev.partition_metrics(truth["families"], clusters["families"])
        if partition != result["partition"]:
            raise ValueError("Stored partition differs from raw families")
        partitions.append(partition)
        for role, collection in (("context", cohorts), ("base", base_cohorts)):
            predictions = load(directory / f"{role}_predictions.json")
            if dev.score(truth, predictions) != result[f"{role}_score"]:
                raise ValueError("Stored score differs from raw rows")
            collection.append({"cohort_id": f"{split}_{config['index']:02d}", "truth": truth["pairs"], "predictions": predictions})
    summaries = {}
    exact = sum(p["exactly_recovered_families"] for p in partitions) / (8 * len(configs))
    missing = sum(len(p["missing_samples"]) for p in partitions)
    for role, collection in (("context", cohorts), ("base", base_cohorts)):
        truth, rows, scores = {"pairs": []}, [], []
        for index, cohort in enumerate(collection):
            scores.append(dev.score({"pairs": cohort["truth"]}, cohort["predictions"]))
            for source, target in ((cohort["truth"], truth["pairs"]), (cohort["predictions"], rows)):
                target.extend({**row, "indID1": f"c{index}_{row['indID1']}", "indID2": f"c{index}_{row['indID2']}"} for row in source)
        overall = dev.score(truth, rows)
        rng, bootstrap = np.random.default_rng(20301120), []
        for _ in range(2000):
            correct, counts = Counter(), Counter()
            for index in rng.integers(0, len(configs), len(configs)):
                correct.update(scores[index]["correct_counts"])
                counts.update(scores[index]["counts"])
            bootstrap.append(np.mean([correct[label] / counts[label] for label in dev.LABELS]))
        overall["cohort_bootstrap_balanced_interval95"] = np.quantile(bootstrap, [.025, .975]).tolist()
        checks = gate_checks(overall, exact, missing, plan["gates"])
        summaries[role] = {"score": overall, "checks": checks, "all_original_gates_passed": all(checks.values())}
        dump(out / f"scoring_{role}.json", {"labels": list(dev.LABELS), "cohorts": collection})
    verify(candidate["bound_inputs"])
    verify(bounds)
    flips = {"base_correct_to_wrong": 0, "base_wrong_to_correct": 0, "strong_base_correct_to_wrong": 0}
    for base, context in zip(base_cohorts, cohorts, strict=True):
        truth = {dev.key(row): row["relation"] for row in context["truth"]}
        baseline = {dev.key(row): row for row in base["predictions"]}
        for row in context["predictions"]:
            previous, label = baseline[dev.key(row)], truth.get(dev.key(row))
            flips["base_correct_to_wrong"] += int(previous["relation"] == label and row["relation"] != label)
            flips["base_wrong_to_correct"] += int(previous["relation"] != label and row["relation"] == label)
            flips["strong_base_correct_to_wrong"] += int(previous["support"] >= .9 and previous["relation"] == label and row["relation"] != label)
    summary = {"split": split, "scope": "Prospective synthetic evidence; no final data or submission authorization",
        **summaries["context"], "base_ablation": summaries["base"], "exact_family_recovery": exact,
        "candidate_sha256": candidate_sha, "plan_sha256": PLAN_SHA, "bound_results": bounds, "cohorts": configs,
        "gates": plan["gates"], "preservation_audit": flips, "output_sha256": {p.name: sha(p) for p in out.iterdir()},
        "submission_ready": False, "independent_root_audit_required": True}
    dump(out / "summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    if not sys.flags.isolated:
        return subprocess.call([sys.executable, "-I", str(Path(__file__).resolve()), *sys.argv[1:]])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("command", choices=("freeze", "run", "summarize"))
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--split", choices=("stress", "confirmation"))
    parser.add_argument("--index", type=int, choices=range(8))
    parser.add_argument("--stress-summary-sha256")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command == "freeze":
        freeze(root)
    elif not args.candidate_sha256 or not args.split:
        parser.error("Frozen candidate hash and registered split required")
    elif args.command == "run":
        if args.index is None:
            parser.error("Cohort index required")
        run(root, args.split, args.index, root / "crossed_context_candidate_v1.json",
            args.candidate_sha256, args.stress_summary_sha256)
    else:
        summarize(root, args.split, root / "crossed_context_candidate_v1.json", args.candidate_sha256)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
