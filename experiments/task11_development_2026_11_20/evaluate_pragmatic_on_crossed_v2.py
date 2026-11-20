"""Fixed old model on all 54 new crossed development cohorts; no fitting."""
import datetime
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def weighted_jaccard(data, predicted_cohorts):
    counts, tp, fp = {}, {}, {}
    def key(row):
        return tuple(sorted((row["indID1"], row["indID2"])))
    for cohort, rows in zip(data, predicted_cohorts, strict=True):
        actual = {key(row): row["relation"] for row in cohort[1]["pairs"]}
        observed = {key(row): row["relation"] for row in rows}
        if len(observed) != len(rows):
            raise ValueError("Duplicate predictions")
        for pair, label in actual.items():
            counts[label] = counts.get(label, 0) + 1
            tp[label] = tp.get(label, 0) + int(observed.get(pair) == label)
        for pair, label in observed.items():
            if actual.get(pair) != label:
                fp[label] = fp.get(label, 0) + 1
    total = sum(counts.values())
    if not total:
        raise ValueError("Empty truth denominator")
    return sum(n * tp.get(label, 0) / (n + fp.get(label, 0)) for label, n in counts.items()) / total


root = Path(sys.argv[1]).resolve()
plan_path = root / "pragmatic_crossed_comparison_plan_v2.json"
plan = load(plan_path)
if len(sys.argv) != 3 or sha(plan_path) != sys.argv[2]:
    raise ValueError("Comparison plan identity mismatch")
for path, digest in plan["bound_inputs"].items():
    if sha(path) != digest:
        raise ValueError(f"Changed comparison binding: {path}")
if plan["bound_inputs"].get(str(Path(__file__).resolve())) != sha(__file__):
    raise ValueError("Unbound comparison script")
if not os.environ.get("SLURM_JOB_ID"):
    raise ValueError("Compute requires Slurm")
if datetime.datetime.now(datetime.UTC) > datetime.datetime.fromisoformat(plan["latest_start_utc"]):
    raise ValueError("Comparison budget deadline passed")
spec = importlib.util.spec_from_file_location("pragmatic_comparison_bindings", root / "pragmatic_context_bindings_v1.py")
bindings = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bindings)
protocol, candidate, ctx, _, infer = bindings.initialize(root, root / "pragmatic_context_candidate_v1.json", plan["candidate_sha256"])
dev = ctx.dev
out = root / "pragmatic_crossed_comparison_v2"
out.mkdir(exist_ok=False)
dev.dump(out / "prospective_plan.json", plan)
dev.dump(out / "job.json", {"job_id": os.environ["SLURM_JOB_ID"]})
crossed = load(root / "crossed_development_plan_v1.json")
configs = [c for c in crossed["cohorts"] if c["split"] == "development"]
if configs != plan["cohorts"] or len(configs) != 54:
    raise ValueError("Changed complete development cohort registry")
data, predictions = [], {"base": [], "context": []}
bounds = {**plan["bound_inputs"], str(plan_path): sha(plan_path)}
map_path = root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map"
for config in configs:
    index = config["index"]
    directory = root / "results_crossed_v1" / f"development_{index:03d}"
    result = load(directory / "result.json")
    if result["cohort"] != config or not result["completed"]:
        raise ValueError("Missing or incorrect development result")
    if load(directory / "protocol.json")["master_sha256"] != sha(root / "crossed_development_plan_v1.json"):
        raise ValueError("Generation plan mismatch")
    dev.require_hashes({str(directory / name): digest for name, digest in result["output_sha256"].items()})
    # Only observed detector features/segments and fixed model enter this call.
    base, context, _, groups = infer(directory, candidate, ctx, map_path)
    dev.dump(out / f"predictions_base_{index:02d}.json", base)
    dev.dump(out / f"predictions_context_{index:02d}.json", context)
    dev.dump(out / f"groups_{index:02d}.json", groups)
    # Truth is consulted only after inference, for unconditional evaluation.
    d = dev.dataset(directory / "result.json")
    data.append(d)
    predictions["base"].append(base)
    predictions["context"].append(context)
    bounds[str(directory / "result.json")] = sha(directory / "result.json")
    print(json.dumps({"completed_cohorts": len(data), "total": 54}), flush=True)
reports = {}
for role, rows in predictions.items():
    score = ctx.aggregate(data, rows)
    reports[role] = {**score, "official_truth_weighted_jaccard": weighted_jaccard(data, rows)}
    dev.dump(out / f"scoring_development_{role}.json", {"labels": list(dev.LABELS), "cohorts": [
        {"cohort_id": f"development_{i:02d}", "truth": d[1]["pairs"], "predictions": predicted}
        for i, (d, predicted) in enumerate(zip(data, rows, strict=True))]})
dev.require_hashes(bounds)
report = {"schema": "task11-fixed-existing-model-crossed-evaluation-v1", "reports": reports,
    "cohorts": 54, "families": 432, "fitted_on_these_cohorts": False, "candidate_sha256": plan["candidate_sha256"],
    "plan_sha256": sha(plan_path), "strict_qualification_passed": False,
    "limitations": plan["limitations"], "inputs_sha256": bounds,
    "outputs_sha256": {p.name: sha(p) for p in sorted(out.iterdir())}, "submission_ready": False}
dev.dump(out / "report.json", report)
print(json.dumps({"completed": True, "reports": reports}))
