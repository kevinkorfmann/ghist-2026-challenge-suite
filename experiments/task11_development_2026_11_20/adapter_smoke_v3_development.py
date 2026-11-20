"""Synthetic-only end-to-end smoke; never use challenge inputs."""
import json
from pathlib import Path
import subprocess
import sys

root = Path(sys.argv[1])
source = root / "source_v5"
sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
from scripts.develop_task11_pedigree import digest  # noqa: E402

candidate = root / "candidate_v5/candidate.json"
project = root.parents[1]
output = root / "adapter_smoke_v3_development"
command = [sys.executable, str(root / "adapter_review_v3/scripts/prepare_task11_candidate.py"),
           "--source", str(source), "--protocol", str(root / "protocol_v5.json"),
           "--protocol-sha256", digest(root / "protocol_v5.json"),
           "--candidate", str(candidate), "--candidate-sha256", digest(candidate),
           "--vcf", str(root / "results_v3/confirmation_00/synthetic.vcf.gz"),
           "--map", str(project / "external/dog-recombination/campbell2016_average.plink.map"),
           "--hapibd", str(project / "external/hap-ibd/hap-ibd.jar"), "--out", str(output),
           "--development-only"]
subprocess.run(command, check=True)


def by_pair(rows):
    return {tuple(sorted((row["indID1"], row["indID2"]))): row for row in rows}


original = by_pair(json.loads((root / "raw_development_v5_v2/predictions_00.json").read_text()))
exported = by_pair(json.loads((output / "predictions.json").read_text()))
if original.keys() != exported.keys():
    raise ValueError("Adapter roster aliasing changed the predicted partition")
changed = sum(original[pair]["relation"] != exported[pair]["relation"] for pair in original)
max_probability_delta = max(abs(original[pair]["probabilities"][label] - exported[pair]["probabilities"][label])
                            for pair in original for label in original[pair]["probabilities"])
report = {"synthetic_input_only": True, "same_biological_cohort_reused": True, "development_only": True,
          "pairs": len(exported), "label_changes": changed,
          "maximum_probability_delta": max_probability_delta,
          "all_passed": changed == 0 and max_probability_delta < 1e-12,
          "adapter_audit_sha256": digest(output / "adapter_audit.json"),
          "reference_predictions_sha256": digest(root / "raw_development_v5_v2/predictions_00.json"),
          "smoke_script_sha256": digest(Path(__file__)), "submission_count": 0}
(root / "adapter_smoke_report_v3_development.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
if not report["all_passed"]:
    raise SystemExit(1)
