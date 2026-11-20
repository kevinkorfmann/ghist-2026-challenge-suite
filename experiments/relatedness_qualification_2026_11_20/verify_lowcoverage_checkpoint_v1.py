"""Verify low-coverage evidence integrity and truth-control score arithmetic."""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import tarfile


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "experiments/lowcoverage_2026_11_20/evidence_checkpoint_v1.tar.gz"


def main():
    digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    assert digest == "cfe3bf20fec866e4d1eac9458630793b863068a386ef1df670990f9ea854abfc"
    with tarfile.open(ARCHIVE) as archive:
        members = archive.getmembers()
        assert all(not Path(m.name).is_absolute() and ".." not in Path(m.name).parts
                   and not m.issym() and not m.islnk() for m in members)
        files = {m.name: m for m in members if m.isfile()}

        def read(name):
            return archive.extractfile(files[name]).read()

        def load(name):
            return json.loads(read(name))

        manifest = load("artifact_manifest.json")
        assert set(files) == {row["path"] for row in manifest["files"]} | {"artifact_manifest.json"}
        for row in manifest["files"]:
            raw = read(row["path"])
            assert len(raw) == row["bytes"]
            assert hashlib.sha256(raw).hexdigest() == row["sha256"], row["path"]
        controls = {}
        for task in ("bottleneck", "split_isolation"):
            estimate = load(f"full_truth_{task}_v5/estimate.json")
            score = load(f"full_truth_{task}_v5/scored.json")
            simulation_path = f"full_{task}/synthetic/simulation.json"
            simulation = load(simulation_path)
            assert hashlib.sha256(read(simulation_path)).hexdigest() == score["simulation_sha256"]
            assert score["truth"] == simulation["truth"]
            assert score["estimate"] == estimate["targets"]
            assert simulation["length"] == estimate["callable_length"] == 10_000_000
            errors = {key: (estimate["targets"][key] - truth) / truth for key, truth in score["truth"].items()}
            assert all(math.isclose(value, score["relative_errors"][key], abs_tol=1e-14)
                       for key, value in errors.items())
            parameters, theta = estimate["best_fit"]["parameters"], estimate["best_fit"]["theta"]
            ancestral = theta / (4 * simulation["mutation_rate"] * estimate["callable_length"])
            assert math.isclose(ancestral, estimate["ancestral_population_size"])
            if task == "bottleneck":
                converted = {"generations": 2 * ancestral * parameters[1], "post_decline_fraction": parameters[0]}
                bounds = [(.001, .9999), (.0001, 10)]
            else:
                converted = {"pop1_population_size": ancestral * parameters[0],
                             "pop2_population_size": ancestral * parameters[1],
                             "generations_since_split": 2 * ancestral * parameters[2]}
                bounds = [(.05, 20), (.05, 20), (.0001, 10)]
            assert all(math.isclose(value, estimate["targets"][key]) for key, value in converted.items())
            runs = estimate["restarts"]
            best = min(runs, key=lambda row: row["negative_log_likelihood"])
            matching = set()
            for row in runs:
                boundary = any(abs(math.log(p / lo)) < 1e-3 or abs(math.log(p / hi)) < 1e-3
                               for p, (lo, hi) in zip(row["parameters"], bounds))
                assert boundary == row["at_boundary"]
                if (row["success"] and not boundary
                        and abs(row["negative_log_likelihood"] - best["negative_log_likelihood"]) < .01
                        and max(abs(math.log(a / b)) for a, b in zip(row["parameters"], best["parameters"])) < .01):
                    matching.add(row["method"])
            assert matching == {"L-BFGS-B", "Powell"} and best["success"] and not best["at_boundary"]
            controls[task] = {"targets": estimate["targets"], "relative_errors": errors,
                              "absolute_scaling_independently_verified": True,
                              "optimizer_agreement_recomputed": True, "restarts": len(runs),
                              "all_targets_within_30percent": all(abs(error) <= .3 for error in errors.values()),
                              "published_demography_root_sum_squared_relative_error": math.sqrt(sum(error**2 for error in errors.values())),
                              "genotype_truth_control_only": True}
    record = {"checked_at_utc": datetime.now(timezone.utc).isoformat(), "archive_sha256": digest,
              "manifest_bound_files_verified": len(manifest["files"]), "total_files_including_manifest": len(files),
              "controls": controls, "release_authorized": False,
              "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Archive integrity, score arithmetic, absolute conversion and optimization-record consistency; no independent raw-BAM fitting rerun or scientific release authorization"}
    output = ROOT / "experiments/relatedness_qualification_2026_11_20/lowcoverage_checkpoint_v1_audit.json"
    with output.open("x") as handle:
        handle.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps(record))


if __name__ == "__main__":
    main()
