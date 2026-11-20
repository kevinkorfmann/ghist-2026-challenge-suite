"""Truth-free CSV adapter for the separately frozen crossed Task11 context model.

Only the VCF roster, observed genotypes, dog map, detector and bound models enter
inference. Qualification summaries are checked separately from the input VCF.
This file never uploads or authorizes final challenge access.
"""
import argparse
import importlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path


def require_qualification(args, runner):
    evidence = {}
    plan = runner.load(args.root.resolve() / "crossed_development_plan_v1.json")
    for split in ("stress", "confirmation"):
        path, expected = getattr(args, f"{split}_summary"), getattr(args, f"{split}_summary_sha256")
        if path is not None:
            runner.verify({str(path): expected})
            value = runner.load(path)
            if (value.get("split") != split or value.get("candidate_sha256") != args.candidate_sha256 or
                    value.get("plan_sha256") != runner.PLAN_SHA or
                    set(value.get("checks", {})) != runner.CHECKS or
                    value.get("all_original_gates_passed") is not True or
                    any(v is not True for v in value["checks"].values()) or
                    value.get("cohorts") != plan["reserved_future_cohorts"][split]):
                raise ValueError(f"Invalid candidate-bound {split} evidence")
            runner.verify(value["bound_results"])
            runner.verify({str(path.parent / name): value for name, value in value["output_sha256"].items()})
            evidence[split] = expected
    passed = set(evidence) == {"stress", "confirmation"}
    if not passed and not args.development_only:
        raise ValueError("Both prospective stress and confirmation must pass; synthetic smoke requires --development-only")
    return {"stress_and_confirmation_gates_passed": passed, "development_only": args.development_only,
            "evidence_sha256": evidence, "scope": "Synthetic evidence only; root independent audit and final-data authorization remain separate"}


def run(args, runner):
    root = args.root.resolve()
    plan, protocol, candidate, ctx, adapter, infer = runner.initialize(
        root, root / "crossed_context_candidate_v1.json", args.candidate_sha256)
    if candidate["bound_inputs"].get(str(Path(__file__).resolve())) != runner.sha(__file__):
        raise ValueError("Context CSV adapter is not bound to the candidate")
    qualification = require_qualification(args, runner)
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Heavy genotype inference requires Slurm")
    project = root.parents[1]
    map_path = project / "external/dog-recombination/campbell2016_average.plink.map"
    hapibd = project / "external/hap-ibd/hap-ibd.jar"
    for path in (map_path, hapibd):
        if str(path) not in protocol["inputs"]:
            raise ValueError("Unbound map/detector")
    dev, np = ctx.dev, ctx.np
    validator = importlib.import_module("ghist_submission.validate")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    input_hash = runner.sha(args.vcf)
    audit = {"schema": "task11-crossed-context-csv-adapter-v1", "candidate_sha256": args.candidate_sha256,
        "plan_sha256": runner.PLAN_SHA, "protocol_sha256": candidate["protocol_sha256"],
        "adapter_sha256": runner.sha(__file__), "input_vcf": str(args.vcf.resolve()),
        "input_vcf_sha256": input_hash, "truth_used": False, "qualification": qualification,
        "submission_ready": False, "submission_count": 0, "completed": False,
        "job_id": os.environ["SLURM_JOB_ID"], "context_parameters": plan["context_parameters"]}
    try:
        roster = adapter.stage_vcf(args.vcf, out / "synthetic.vcf.gz")
        runner.dump(out / "roster.json", roster)
        inference = dev.detect_cluster_extract(out, map_path, hapibd, roster["internal_samples"], "average")
        for role, count in (("base", 149), ("context", 185)):
            model = dev.joblib.load(candidate[f"{role}_model"])
            if set(model.classes_) != set(protocol["labels"]) or model.n_features_in_ != count:
                raise ValueError(f"Unexpected {role} model labels/features")
        base, rows, extra, grouping = infer(out, candidate, ctx, map_path)
        families = inference["clusters"]["families"]
        rows, families = adapter.restore_predictions(rows, families, roster["internal_to_original"])
        base, _ = adapter.restore_predictions(base, inference["clusters"]["families"], roster["internal_to_original"])
        expected = validator.within_family_pairs(families, roster["samples"])
        adapter.export_csv(out / "task11_candidate.csv", rows)
        validation = validator.validate_submission(out / "task11_candidate.csv", 11, args.dataset, expected)
        runner.dump(out / "predictions.json", rows)
        runner.dump(out / "base_predictions.json", base)
        runner.dump(out / "families.json", {"families": families, "supplied_family_count": 8})
        runner.dump(out / "validation.json", validation)
        runner.dump(out / "grouping_audit.json", grouping)
        np.savez_compressed(out / "context_features.npz", features=extra)
        runner.verify({**candidate["bound_inputs"], str(args.vcf): input_hash})
        adapter.verify_source(root / "source_v5", protocol["source_sha256"])
        dev.require_runtime(protocol)
        require_qualification(args, runner)
        audit.update(completed=True, samples=len(roster["samples"]), pairs=len(rows), validation=validation,
            stages=inference["stages"], imported_sources=dev.imported_source_tree(),
            base_model_sha256=runner.sha(candidate["base_model"]), context_model_sha256=runner.sha(candidate["context_model"]))
    except Exception as error:
        audit["failure"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        audit["output_sha256"] = {str(p.relative_to(out)): runner.sha(p) for p in sorted(out.rglob("*")) if p.is_file() and p.name != "adapter_audit.json"}
        runner.dump(out / "adapter_audit.json", audit)
    return audit


def main():
    if not sys.flags.isolated:
        return subprocess.call([sys.executable, "-I", str(Path(__file__).resolve()), *sys.argv[1:]])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    for name in ("vcf", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("candidate-sha256", "runner-sha256"):
        parser.add_argument("--" + name, required=True)
    for split in ("stress", "confirmation"):
        parser.add_argument(f"--{split}-summary", type=Path)
        parser.add_argument(f"--{split}-summary-sha256")
    parser.add_argument("--development-only", action="store_true")
    parser.add_argument("--dataset", choices=("testing", "final"), default="testing")
    args = parser.parse_args()
    if args.development_only and args.dataset != "testing":
        parser.error("Development-only smoke requires testing format context")
    import hashlib
    runner_path = args.root.resolve() / "run_crossed_future_v1.py"
    if hashlib.sha256(runner_path.read_bytes()).hexdigest() != args.runner_sha256:
        raise ValueError("Reserved-cohort runner hash mismatch before import")
    spec = importlib.util.spec_from_file_location("crossed_inference_bindings", runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    audit = run(args, runner)
    print(json.dumps({"completed": audit["completed"], "pairs": audit["pairs"], "qualification": audit["qualification"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
