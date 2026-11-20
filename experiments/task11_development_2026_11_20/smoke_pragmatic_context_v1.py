"""Retained synthetic VCF smoke with simultaneous sample aliases/order changes.

One complete detector-to-CSV rerun, not additional biological validation or
160 classification reruns. Direct fit predictions are consulted only afterward.
"""
import gzip
import hashlib
import importlib.util
import json
import os
import random
import subprocess
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def key(row):
    return tuple(sorted((row["indID1"], row["indID2"])))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


root, candidate_sha = Path(sys.argv[1]).resolve(), sys.argv[2]
if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Synthetic detector smoke requires Slurm")
candidate_path = root / "pragmatic_context_candidate_v1.json"
if sha(candidate_path) != candidate_sha:
    raise ValueError("Candidate hash mismatch")
candidate = load(candidate_path)
out = root / "pragmatic_context_smoke_v1"
out.mkdir(exist_ok=False)
vcf = root / "results_v5/calibration_00/synthetic.vcf.gz"
original_hash = sha(vcf)
permuted = out / "permuted_synthetic.vcf.gz"
randomness = random.Random(20301122)
mapping, order, records = {}, None, 0
with gzip.open(vcf, "rt") as reader, gzip.open(permuted, "wt") as writer:
    for line in reader:
        if line.startswith("##"):
            writer.write(line)
        elif line.startswith("#CHROM\t"):
            fields = line.rstrip("\n").split("\t")
            samples = fields[9:]
            order = list(range(len(samples)))
            randomness.shuffle(order)
            aliases = [f"synthetic_alias_{randomness.getrandbits(64):016x}" for _ in samples]
            if len(set(aliases)) != len(aliases) or order == list(range(len(samples))):
                raise ValueError("Smoke did not create unique aliases and a nontrivial permutation")
            mapping = dict(zip(aliases, samples, strict=True))
            writer.write("\t".join(fields[:9] + [aliases[i] for i in order]) + "\n")
        else:
            fields = line.rstrip("\n").split("\t")
            if order is None or len(fields) != 9 + len(order):
                raise ValueError("Invalid synthetic VCF layout")
            writer.write("\t".join(fields[:9] + [fields[9 + i] for i in order]) + "\n")
            records += 1
write(out / "alias_and_order_map.json", {"alias_to_original": mapping, "column_order": order,
    "seed": 20301122, "variant_records": records, "original_vcf_sha256": original_hash})
adapter = root / "prepare_pragmatic_context_v1.py"
runner = root / "pragmatic_context_bindings_v1.py"
if candidate["bound_inputs"].get(str(adapter)) != sha(adapter) or candidate["bound_inputs"].get(str(runner)) != sha(runner):
    raise ValueError("Unbound adapter/runner")
command = [sys.executable, "-I", str(adapter), str(root), "--vcf", str(permuted), "--out", str(out / "adapter"),
           "--candidate-sha256", candidate_sha, "--runner-sha256", sha(runner), "--dataset", "testing"]
with (out / "adapter.stdout").open("w") as stdout, (out / "adapter.stderr").open("w") as stderr:
    subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
spec = importlib.util.spec_from_file_location("pragmatic_smoke_bindings", runner)
bindings = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bindings)
_, _, ctx, _, infer = bindings.initialize(root, candidate_path, candidate_sha)
original = root / "results_v5/calibration_00"
raw = load(original / "result.json")
bindings.verify({str(original / name): value for name, value in raw["output_sha256"].items()})
direct_base, direct_context, _, _ = infer(original, candidate, ctx,
    root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map")
fit = root / "context_diagnostic_v1"
report = load(fit / "report.json")
comparison, bounds = {}, {}
direct = {"base": direct_base, "context": direct_context}
write(out / "direct_reference_predictions.json", direct)
for role, filename in (("base", "base_predictions.json"), ("context", "predictions.json")):
    reference = {"predictions": direct[role]}
    reference_path = out / "direct_reference_predictions.json"
    if role == "context":
        historic_path = fit / "scoring_calibration.json"
        if report["outputs_sha256"][historic_path.name] != sha(historic_path):
            raise ValueError("Historic context reference changed")
        historic = load(historic_path)["cohorts"][0]
        if historic["cohort_id"] != "calibration_00":
            raise ValueError("Unexpected historic context cohort")
        by_pair = {key(row): row for row in historic["predictions"]}
        if len(by_pair) != len(direct_context) or set(by_pair) != {key(row) for row in direct_context}:
            raise ValueError("Direct context changed historical pair coverage")
        for row in direct_context:
            old = by_pair[key(row)]
            if old["relation"] != row["relation"] or any(abs(old["probabilities"][label] - probability) > 1e-12
                    for label, probability in row["probabilities"].items()):
                raise ValueError("Direct context no longer reproduces historical inference")
        bounds[str(historic_path)] = sha(historic_path)
    actual = load(out / "adapter" / filename)
    actual = [{**row, "indID1": mapping[row["indID1"]], "indID2": mapping[row["indID2"]]} for row in actual]
    observed, expected = {key(r): r for r in actual}, {key(r): r for r in reference["predictions"]}
    if len(observed) != len(actual) or len(expected) != len(reference["predictions"]) or set(observed) != set(expected):
        raise ValueError("Prediction pair set changed under alias/order permutation")
    changes = sum(observed[pair]["relation"] != expected[pair]["relation"] for pair in expected)
    delta = max(abs(observed[pair]["probabilities"][label] - expected[pair]["probabilities"][label])
                for pair in expected for label in expected[pair]["probabilities"])
    comparison[role] = {"pairs": len(expected), "label_changes": changes, "maximum_probability_delta": delta,
                        "passed": changes == 0 and delta <= 1e-12}
    if not comparison[role]["passed"]:
        raise ValueError(f"{role} alias/order parity failed")
    bounds[str(reference_path)] = sha(reference_path)
clusters = load(root / "results_v5/calibration_00/clusters.json")["families"]
restored_families = [[mapping[s] for s in family] for family in load(out / "adapter/families.json")["families"]]
if {frozenset(f) for f in clusters} != {frozenset(f) for f in restored_families}:
    raise ValueError("Cluster partition changed under alias/order permutation")
if sha(vcf) != original_hash:
    raise ValueError("Synthetic source changed")
audit = {"scope": "One retained synthetic cohort with simultaneous alias and column-order permutation; operational parity only",
    "candidate_sha256": candidate_sha, "protocol_sha256": candidate["protocol_sha256"],
    "source_vcf_sha256": original_hash, "permuted_vcf_sha256": sha(permuted),
    "comparisons": comparison, "family_partition_unchanged": True, "full_classification_reruns": 1,
    "adapter_audit_sha256": sha(out / "adapter/adapter_audit.json"), "reference_bounds": bounds,
    "smoke_script_sha256": sha(__file__), "job_id": os.environ["SLURM_JOB_ID"],
    "submission_ready": False, "new_biological_validation": False,
    "output_sha256": {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()}}
write(out / "smoke_audit.json", audit)
print(json.dumps(audit, indent=2))
