"""Retained synthetic map-crop identity and full unchanged-adapter parity."""
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
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


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def key(row):
    return tuple(sorted((row["indID1"], row["indID2"])))


root, plan_sha = Path(sys.argv[1]).resolve(), sys.argv[2]
plan_path = root / "map_bounds_staging_plan_v1.json"
if sha(plan_path) != plan_sha or not os.environ.get("SLURM_JOB_ID"):
    raise ValueError("Expected bound plan and Slurm")
plan = load(plan_path)
for path, digest in plan["bound_inputs"].items():
    if sha(path) != digest:
        raise ValueError(f"Changed smoke input: {path}")
spec = importlib.util.spec_from_file_location("map_bounds_stage", root / "stage_task11_map_bounds_v1.py")
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)
out = root / "map_bounds_smoke_v1"
out.mkdir(exist_ok=False)
vcf = root / "results_v5/calibration_00/synthetic.vcf.gz"
map_path = root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map"
original_result = load(root / "results_v5/calibration_00/result.json")
if sha(vcf) != original_result["output_sha256"]["synthetic.vcf.gz"]:
    raise ValueError("Changed retained synthetic VCF")
identity = stage.stage(vcf, map_path, out / "identity.vcf.gz", plan["map_sha256"])
if identity["excluded_records"] != 0 or identity["input_variant_lines_sha256"] != identity["retained_variant_lines_sha256"]:
    raise ValueError("Retained synthetic input is not an identity under map crop")
write(out / "identity_stage_audit.json", identity)
# Add one copied genotype row just above each map endpoint, in chromosome order.
bounds = stage.bounds_from_map(map_path)
last, added = None, 0
spiked = out / "synthetic_with_out_of_map_markers.vcf.gz"
with gzip.open(vcf, "rb") as reader, gzip.open(spiked, "wb") as writer:
    for line in reader:
        if line.startswith(b"#"):
            writer.write(line)
            continue
        fields = line.rstrip(b"\r\n").split(b"\t")
        if last is not None and fields[0] != last[0]:
            extra = list(last)
            extra[1] = str(bounds[last[0].decode()][1] + 1).encode()
            writer.write(b"\t".join(extra) + b"\n")
            added += 1
        writer.write(line)
        last = fields
    extra = list(last)
    extra[1] = str(bounds[last[0].decode()][1] + 1).encode()
    writer.write(b"\t".join(extra) + b"\n")
    added += 1
cropped = out / "cropped_spiked.vcf.gz"
crop = stage.stage(spiked, map_path, cropped, plan["map_sha256"])
if added != 38 or crop["excluded_records"] != added or crop["output_vcf_sha256"] != identity["output_vcf_sha256"]:
    raise ValueError("Out-of-map spike crop did not reproduce the exact identity output")
write(out / "spiked_stage_audit.json", crop)
command = [sys.executable, "-I", str(root / "prepare_pragmatic_context_v1.py"), str(root),
    "--vcf", str(cropped), "--out", str(out / "adapter"), "--dataset", "testing",
    "--candidate-sha256", plan["candidate_sha256"], "--runner-sha256", plan["binding_script_sha256"]]
with (out / "adapter.stdout").open("w") as stdout, (out / "adapter.stderr").open("w") as stderr:
    subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
direct = load(root / "pragmatic_context_smoke_v1/direct_reference_predictions.json")
comparison = {}
for role, filename in (("base", "base_predictions.json"), ("context", "predictions.json")):
    expected, actual_rows = {key(row): row for row in direct[role]}, load(out / "adapter" / filename)
    actual = {key(row): row for row in actual_rows}
    if set(actual) != set(expected) or len(actual) != len(actual_rows):
        raise ValueError("Map staging changed pair coverage")
    changed = sum(actual[p]["relation"] != expected[p]["relation"] for p in expected)
    delta = max(abs(actual[p]["probabilities"][label] - value) for p in expected for label, value in expected[p]["probabilities"].items())
    comparison[role] = {"pairs": len(expected), "label_changes": changed, "maximum_probability_delta": delta,
                        "passed": changed == 0 and delta <= 1e-12}
    if not comparison[role]["passed"]:
        raise ValueError("Map staging changed fixed-model predictions")
expected_families = load(root / "results_v5/calibration_00/clusters.json")["families"]
actual_families = load(out / "adapter/families.json")["families"]
if {frozenset(f) for f in expected_families} != {frozenset(f) for f in actual_families}:
    raise ValueError("Map staging changed family partition")
for path, digest in plan["bound_inputs"].items():
    if sha(path) != digest:
        raise ValueError("Bound source changed during smoke")
audit = {"schema": "task11-map-bound-operational-smoke-v1", "plan_sha256": plan_sha,
    "candidate_sha256": plan["candidate_sha256"], "identity_exclusions": identity["excluded_records"],
    "added_out_of_map_records": added, "removed_out_of_map_records": crop["excluded_records"],
    "cropped_spike_exactly_matches_identity_gzip": True, "retained_genotypes_and_headers_unchanged": True,
    "comparisons": comparison, "family_partition_unchanged": True, "full_classification_reruns": 1,
    "new_biological_validation": False, "job_id": os.environ["SLURM_JOB_ID"], "submission_ready": False,
    "output_sha256": {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()}}
write(out / "smoke_audit.json", audit)
print(json.dumps(audit, indent=2))
