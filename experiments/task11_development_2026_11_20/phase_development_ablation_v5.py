"""Matched phase/no-phase development ablation; exposed cohorts are not confirmation."""
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
source = root / "source_v5"
sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
from scripts import develop_task11_pedigree as dev  # noqa: E402

protocol = dev.load(root / "protocol_v5.json")
dev.require_runtime(protocol)
dev.require_hashes({str(source / path): h for path, h in protocol["source_sha256"].items()})
candidate = dev.load(root / "candidate_v5/candidate.json")
baseline = dev.load(root / "candidate_v3/candidate.json")
dev.require_hashes(candidate["bound_inputs"])
dev.require_hashes(baseline["bound_inputs"])
phase_model = dev.joblib.load(candidate["model_path"])
base_model = dev.joblib.load(baseline["model_path"])
matched = []
for split in ("train", "calibration"):
    for index in range(dev.COUNTS[split]):
        old = dev.dataset(root / "results_v3" / f"{split}_{index:02d}" / "result.json")
        new = dev.dataset(root / "results_v5" / f"{split}_{index:02d}" / "result.json")
        if old[1] != new[1] or old[2] != new[2]:
            raise ValueError("Ablation does not have matched whole-cohort truth/features")
        dev.np.testing.assert_allclose(old[3], new[3][:, :-18], rtol=1e-12, atol=1e-12)
        matched.append({"split": split, "index": index, "base_features_equal": True})

reports, truth_all, base_all, phase_all = [], {"pairs": []}, [], []
phase_output = root / "phase_ablation_v5"
phase_output.mkdir(exist_ok=False)
for index in range(8):
    cohort = root / "results_v3" / f"confirmation_{index:02d}"
    result, truth, rows, base_features, *_ = dev.dataset(cohort / "result.json")
    phase_features = dev.haplotype_features(cohort / "segments.ibd.gz", cohort / "segments.hbd.gz",
                    root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map",
                    [(r["indID1"], r["indID2"]) for r in rows])
    with_phase = dev.np.concatenate([base_features, phase_features], axis=1)
    phase_path = phase_output / f"exposed_cohort_{index:02d}.npz"
    dev.np.savez_compressed(phase_path, features=phase_features)
    base_rows = dev.predict(rows, base_features, base_model, baseline["temperature"], baseline["sibling_pooling"])
    phase_rows = dev.predict(rows, with_phase, phase_model, candidate["temperature"], candidate["sibling_pooling"])
    reports.append({"cohort": result["cohort"], "base": dev.score(truth, base_rows),
                    "phase": dev.score(truth, phase_rows),
                    "phase_features_path": str(phase_path), "phase_features_sha256": dev.digest(phase_path),
                    "zero_endpoint_denominators": {str(t): int((phase_features[:, i*6+2:i*6+4] == 0).sum()) for i, t in enumerate(dev.HBD)},
                    "retained_fraction_quantiles": {str(t): dev.np.quantile(phase_features[:, i*6+4:i*6+6], [0,.25,.5,.75,1]).tolist() for i,t in enumerate(dev.HBD)}})
    for original, target in ((truth["pairs"], truth_all["pairs"]), (base_rows, base_all), (phase_rows, phase_all)):
        target.extend({**row, "indID1": f"c{index}_{row['indID1']}", "indID2": f"c{index}_{row['indID2']}"} for row in original)
report = {"scope": "Development-only reuse of eight previously exposed failed-confirmation cohorts; neither fitted nor used as fresh confirmation",
          "matched_base_features": matched, "base": dev.score(truth_all, base_all), "phase": dev.score(truth_all, phase_all),
          "cohorts": reports, "protocol_sha256": dev.digest(root / "protocol_v5.json"),
          "candidate_sha256": dev.digest(root / "candidate_v5/candidate.json"),
          "baseline_candidate_sha256": dev.digest(root / "candidate_v3/candidate.json"),
          "script_sha256": dev.digest(Path(__file__)), "submission_ready": False}
dev.dump(root / "phase_development_ablation_v5.json", report)
print(json.dumps({"base": report["base"], "phase": report["phase"]}))
