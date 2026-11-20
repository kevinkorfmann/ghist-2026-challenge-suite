"""Independently recompute full-width BAM-development gates from retained evidence."""
import hashlib
import json
import math
import re
import tarfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/lowcoverage_2026_11_20"
ARCHIVE = BASE / "evidence_checkpoint_v2.tar.gz"
ARCHIVE_SHA = "143baed0746831ded2b90abe7d65abffdca63d8eb81e1eb51bbdd6a49e6b11f0"
POLICY_SHA = "6bb7811d87d0028dfa5d817886e5d26dbe171585e7939b4fc027e8e28d4415f6"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/lowcoverage_checkpoint_v2_audit.json"


def close(a, b):
    assert math.isclose(a, b, abs_tol=1e-12, rel_tol=1e-12), (a, b)


def optimizer_agreement(estimate, task):
    bounds = [(.001, .9999), (.0001, 10)] if task == "bottleneck" else [(.05, 20), (.05, 20), (.0001, 10)]
    runs = estimate["restarts"]
    assert len(runs) == 16
    assert Counter(row["method"] for row in runs) == {"L-BFGS-B": 8, "Powell": 8}
    best = min(runs, key=lambda row: row["negative_log_likelihood"])
    assert best == estimate["best_fit"]
    methods = set()
    for row in runs:
        assert math.isfinite(row["negative_log_likelihood"])
        boundary = any(abs(math.log(p / lo)) < 1e-3 or abs(math.log(p / hi)) < 1e-3
                       for p, (lo, hi) in zip(row["parameters"], bounds, strict=True))
        assert boundary == row["at_boundary"]
        if (row["success"] and not boundary
                and abs(row["negative_log_likelihood"] - best["negative_log_likelihood"]) < .01
                and max(abs(math.log(a / b)) for a, b in zip(row["parameters"], best["parameters"], strict=True)) < .01):
            methods.add(row["method"])
    passed = best["success"] and not best["at_boundary"] and methods == {"L-BFGS-B", "Powell"}
    assert passed == estimate["optimizer_agreement"]["passed"]
    return passed


def main():
    assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() == ARCHIVE_SHA
    with tarfile.open(ARCHIVE) as archive:
        members = archive.getmembers()
        assert len({m.name for m in members}) == len(members)
        assert all(not Path(m.name).is_absolute() and ".." not in Path(m.name).parts
                   and not m.issym() and not m.islnk() for m in members)
        files = {m.name: m for m in members if m.isfile()}

        def read(name):
            return archive.extractfile(files[name]).read()

        def load(name):
            return json.loads(read(name))

        def digest(name):
            return hashlib.sha256(read(name)).hexdigest()

        manifest = load("artifact_manifest.json")
        assert set(files) == {row["path"] for row in manifest["files"]} | {"artifact_manifest.json"}
        for row in manifest["files"]:
            assert len(read(row["path"])) == row["bytes"]
            assert digest(row["path"]) == row["sha256"], row["path"]
        assert digest("winsfs_full_development_v1.json") == POLICY_SHA
        source = load("source_v7/source_manifest.json")
        for name, expected in source.items():
            assert digest("source_v7/" + name) == expected
        claimed = load("full_development_gate_audit_v7.json")
        tasks = {}
        for task in ("bottleneck", "split_isolation"):
            prefix = f"window_full_{task}_v7"
            simulation_path = f"full_{task}/synthetic/simulation.json"
            simulation = load(simulation_path)
            scored = load(prefix + "/scored.json")
            binding = load(prefix + "/execution_binding.json")
            initial_path = f"full_execution_initial_{task}_v7.json"
            initial = load(initial_path)
            assert {key: binding[key] for key in initial} == initial
            assert binding["source_and_input_pre_post_verification_passed"]
            for field, path in (
                ("initial_record_sha256", initial_path),
                ("scored_result_sha256", prefix + "/scored.json"),
                ("likelihood_protocol_sha256", prefix + "/safs/protocol.json"),
                ("window_protocol_sha256", prefix + "/window_comparison/protocol.json"),
                ("source_manifest_sha256", "source_v7/source_manifest.json"),
                ("launch_script_sha256", "winsfs_full_v7.sbatch"),
                ("environment_record_sha256", "environment_v7.json"),
            ):
                assert binding[field] == digest(path), field
            assert binding["source_hashes_verified"] == source
            assert binding["protocol_sha256"] == POLICY_SHA
            assert binding["input_audit"] == load(prefix + "/safs/inputs.json")
            assert scored["simulation_sha256"] == digest(simulation_path)
            assert scored["truth"] == simulation["truth"]
            assert not scored["testing_data_access"] and not scored["final_data_access"]
            likelihood = load(prefix + "/safs/protocol.json")
            window = load(prefix + "/window_comparison/protocol.json")
            assert likelihood["angsd_sha256"] == binding["external_tool_sha256"]["external/angsd/angsd"]
            assert window["winsfs_sha256"] == binding["external_tool_sha256"]["external/winsfs-release-v0.7.0/winsfs"]
            assert (window["max_epochs"], window["tolerance"], window["blocks"], window["window_size"]) == (50, 1e-4, 500, 100)
            estimates, spectra, errors, epochs, stops, agreement = {}, {}, {}, {}, {}, {}
            for label in ("seed1", "seed2", "continued"):
                folder = prefix + "/window_comparison/"
                estimate = load(folder + label + ".estimate.json")
                text = read(folder + label + ".folded.sfs").decode()
                header = re.match(r"^#SHAPE=<([0-9/, ]+)>\s*\n", text)
                assert header
                shape = [int(value) for value in re.split(r"[,/ ]+", header[1])]
                assert shape == [n + 1 for n in window["sample_sizes"]]
                values = [float(value) for value in text[header.end():].split()]
                assert len(values) == math.prod(shape) and all(math.isfinite(v) and v >= 0 for v in values)
                exposure = math.fsum(values)
                close(exposure, estimate["callable_length"])
                close(exposure / simulation["length"], estimate["callable_fraction"])
                parameters, theta = estimate["best_fit"]["parameters"], estimate["best_fit"]["theta"]
                ancestral = theta / (4 * simulation["mutation_rate"] * exposure)
                close(ancestral, estimate["ancestral_population_size"])
                targets = ({"generations": 2 * ancestral * parameters[1], "post_decline_fraction": parameters[0]}
                           if task == "bottleneck" else {"pop1_population_size": ancestral * parameters[0],
                           "pop2_population_size": ancestral * parameters[1], "generations_since_split": 2 * ancestral * parameters[2]})
                for key, value in targets.items():
                    close(value, estimate["targets"][key])
                    close(value, scored["comparison"]["estimates"][label][key])
                errors[label] = {key: value / simulation["truth"][key] - 1 for key, value in targets.items()}
                for key, value in errors[label].items():
                    close(value, scored["relative_errors"][label][key])
                log = read(folder + label + ".log").decode()
                epochs[label] = [int(value) for value in re.findall(r"Finished epoch ([0-9]+)", log)]
                stop_lines = re.findall(r"Δ=([-+0-9.e]+) ≤ ([0-9.e-]+)", log)
                stops[label] = bool(epochs[label] and epochs[label][-1] < 50 and stop_lines
                                    and abs(float(stop_lines[-1][0])) <= float(stop_lines[-1][1]) <= 1e-4)
                assert epochs[label] == estimate["winsfs_epochs"]
                agreement[label] = optimizer_agreement(estimate, task)
                estimates[label], spectra[label] = estimate, values
            comparisons = {}
            for label in ("seed2", "continued"):
                a, b = spectra["seed1"], spectra[label]
                comparisons[label] = {
                    "frequency_L1": math.fsum(abs(x / math.fsum(a) - y / math.fsum(b)) for x, y in zip(a, b, strict=True)),
                    "segregating_relative_L1": math.fsum(abs(x - y) for x, y in zip(a[1:], b[1:], strict=True)) / math.fsum(a[1:]),
                    "target_relative_change": {key: estimates[label]["targets"][key] / value - 1
                                               for key, value in estimates["seed1"]["targets"].items()},
                }
                for metric in ("frequency_L1", "segregating_relative_L1"):
                    close(comparisons[label][metric], scored["comparison"]["comparisons"][label][metric])
                for key, value in comparisons[label]["target_relative_change"].items():
                    close(value, scored["comparison"]["comparisons"][label]["target_relative_change"][key])
            gates = {
                "exposure": all(e["callable_fraction"] >= .99 for e in estimates.values()),
                "observed_stopping": all(stops.values()), "two_optimizer_agreement": all(agreement.values()),
                "target_stability": all(abs(x) < .01 for c in comparisons.values() for x in c["target_relative_change"].values()),
                "frequency_L1": all(c["frequency_L1"] < .001 for c in comparisons.values()),
                "segregating_relative_L1": all(c["segregating_relative_L1"] < .05 for c in comparisons.values()),
                "parameter_recovery": all(abs(x) <= .30 for row in errors.values() for x in row.values()),
                "pre_post_binding": True,
            }
            assert gates == claimed["tasks"][task]["gates"]
            tasks[task] = {"gates": gates, "full_bam_development_gates_passed": all(gates.values()),
                           "targets": {label: e["targets"] for label, e in estimates.items()},
                           "relative_errors": errors, "comparisons": comparisons,
                           "optimizer_agreement": agreement, "last_epoch": {label: values[-1] for label, values in epochs.items()},
                           "documented_demography_root_sum_squared_relative_error": {
                               label: math.sqrt(math.fsum(x * x for x in row.values())) for label, row in errors.items()}}
    report = {"verified_at_utc": datetime.now(UTC).isoformat(), "archive_sha256": ARCHIVE_SHA,
              "manifest_bound_files_verified": len(manifest["files"]), "total_files_including_manifest": len(files),
              "tasks": tasks, "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Independent archive/source/provenance checks, raw-spectrum and log diagnostics, optimizer-record and target arithmetic. No independent raw-BAM likelihood rerun; one development cohort per task.",
              "fresh_confirmation_verified": False, "release_authorized": False}
    with OUT.open("x") as handle:
        handle.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({task: {"gates": row["gates"], "targets": row["targets"]["seed1"]} for task, row in tasks.items()}))


if __name__ == "__main__":
    main()
