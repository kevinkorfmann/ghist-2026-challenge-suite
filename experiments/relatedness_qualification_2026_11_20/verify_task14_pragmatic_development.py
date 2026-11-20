"""Independently audit both frozen cap500 Task 14 pragmatic development cases."""
import argparse
import hashlib
import json
import math
import re
import statistics
import tarfile
from datetime import UTC, datetime
from pathlib import Path

from verify_lowcoverage_checkpoint_v2 import close, optimizer_agreement


FREEZE_SHA = "7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb"
TARGETS = ("pop1_population_size", "pop2_population_size", "generations_since_split")


def audit(args):
    assert hashlib.sha256(args.archive.read_bytes()).hexdigest() == args.sha256
    with tarfile.open(args.archive) as archive:
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
        assert digest("task14_pragmatic_freeze_v1.json") == FREEZE_SHA
        frozen = load("task14_pragmatic_freeze_v1.json")
        for name, expected in frozen["source_files"].items():
            assert digest("source_v12/" + name) == expected
        externally_bound = {}
        for name, expected in frozen["bindings"].items():
            if name in files:
                assert digest(name) == expected, name
            else:
                assert name.startswith("external/"), name
                externally_bound[name] = expected
        original = load("task14_confirmation_freeze_v1.json")
        assert len(original["confirmation_designs"]) == 3
        cohorts, primary, noise_errors, seeds = {}, [], {0: [], 1: []}, []
        for regime, design in enumerate(original["confirmation_designs"][:1]):
            for replicate in range(2):
                stem = f"split_r{regime}_p{replicate}_v1"
                case = "full_confirmation_" + stem
                output = "window_pragmatic_development_" + stem
                simulation_name = case + "/synthetic/simulation.json"
                simulation = load(simulation_name)
                expected_design = dict(design)
                seed = expected_design.pop("seeds")[replicate]
                expected_design["seed"] = seed
                seeds.append(seed)
                assert simulation["design"] == expected_design
                assert simulation["length"] == 10_000_000
                close(simulation["mutation_rate"], 1e-7)
                noise = original["noise_profiles"][f"replicate{replicate}"]
                assert simulation["quality"] == noise["quality"]
                assert simulation["individual_depth_lognormal_sigma"] == noise["depth_sigma"]
                assert not simulation["testing_data_access"] and not simulation["final_data_access"]
                truth = dict(zip(TARGETS, (design["size1"], design["size2"], design["generations"]), strict=True))
                for key, value in truth.items():
                    close(value, simulation["truth"][key])
                initial = load(output + "/initial_binding.json")
                execution = load(output + "/execution_binding.json")
                for record in (initial, execution):
                    assert record["freeze_sha256"] == FREEZE_SHA
                    assert record["runtime_source_tools_formatter_verified"] is True
                    assert (record["regime"], record["replicate"], record["mode"]) == (regime, replicate, "development")
                    assert record["wrapper_sha256"] == frozen["source_files"]["scripts/lowcoverage_pragmatic.py"]
                clock = initial["clock_guard"]
                assert clock["stage"] == "development"
                assert datetime.fromisoformat(clock["actual_start_utc"]) <= datetime.fromisoformat("2026-11-20T20:45:00+00:00")
                assert initial["clock_guard"] == execution["clock_guard"]
                assert execution["initial_binding_sha256"] == digest(output + "/initial_binding.json")
                assert execution["source_runtime_input_pre_post_verified"] is True
                assert execution["synthetic_formatter_schema_passed"] is True
                binding = initial["input_binding"]
                assert execution["input_binding"] == binding
                previous = "window_confirmation_" + stem
                assert load(previous + "/execution_binding.json")["input_binding"] == binding
                assert load(case + "/prepare_binding.json")["input_binding"] == binding
                assert binding["simulation_sha256"] == digest(simulation_name)
                assert binding["audit"] == load(output + "/safs/inputs.json")
                for name, expected in execution["output_bindings"].items():
                    assert digest(output + "/" + name) == expected
                old_safs = load(previous + "/window_comparison/protocol.json")["saf_files"]
                new_safs = load(output + "/window_comparison/protocol.json")["saf_files"]
                assert old_safs == new_safs
                scored = load(output + "/scored.json")
                assert scored["simulation_sha256"] == digest(simulation_name)
                assert scored["truth"] == simulation["truth"]
                assert not scored["testing_data_access"] and not scored["final_data_access"]
                likelihood = load(output + "/safs/protocol.json")
                window = load(output + "/window_comparison/protocol.json")
                assert likelihood["angsd_sha256"] == frozen["bindings"]["external/angsd/angsd"]
                assert window["winsfs_sha256"] == frozen["bindings"]["external/winsfs-release-v0.7.0/winsfs"]
                assert (window["max_epochs"], window["tolerance"], window["blocks"], window["window_size"]) == (500, 1e-6, 500, 500)
                estimates, spectra, errors, stops, agreement, last_epochs = {}, {}, {}, {}, {}, {}
                for label in ("seed1", "seed2", "continued"):
                    prefix = output + "/window_comparison/" + label
                    estimate = load(prefix + ".estimate.json")
                    text = read(prefix + ".folded.sfs").decode()
                    header = re.match(r"^#SHAPE=<([0-9/, ]+)>\s*\n", text)
                    assert header
                    shape = [int(v) for v in re.split(r"[,/ ]+", header[1])]
                    assert shape == [n + 1 for n in window["sample_sizes"]] == [21, 21]
                    values = [float(v) for v in text[header.end():].split()]
                    assert len(values) == math.prod(shape)
                    assert all(math.isfinite(v) and v >= 0 for v in values)
                    exposure = math.fsum(values)
                    close(exposure, estimate["callable_length"])
                    close(exposure / simulation["length"], estimate["callable_fraction"])
                    ancestral = estimate["best_fit"]["theta"] / (4 * simulation["mutation_rate"] * exposure)
                    close(ancestral, estimate["ancestral_population_size"])
                    pars = estimate["best_fit"]["parameters"]
                    targets = dict(zip(TARGETS, (ancestral * pars[0], ancestral * pars[1], 2 * ancestral * pars[2]), strict=True))
                    errors[label] = {key: value / truth[key] - 1 for key, value in targets.items()}
                    for key, value in targets.items():
                        close(value, estimate["targets"][key])
                        close(value, scored["comparison"]["estimates"][label][key])
                        close(errors[label][key], scored["relative_errors"][label][key])
                    log = read(prefix + ".log").decode()
                    epochs = [int(v) for v in re.findall(r"Finished epoch ([0-9]+)", log)]
                    stopping = re.findall(r"Δ=([-+0-9.e]+) ≤ ([0-9.e-]+)", log)
                    assert epochs == estimate["winsfs_epochs"]
                    last_epochs[label] = epochs[-1] if epochs else None
                    stops[label] = bool(epochs and epochs[-1] < 500 and stopping
                                        and abs(float(stopping[-1][0])) <= float(stopping[-1][1]) <= 1e-6)
                    agreement[label] = optimizer_agreement(estimate, "split_isolation")
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
                    for key in ("frequency_L1", "segregating_relative_L1"):
                        close(comparisons[label][key], scored["comparison"]["comparisons"][label][key])
                    for key, value in comparisons[label]["target_relative_change"].items():
                        close(value, scored["comparison"]["comparisons"][label]["target_relative_change"][key])
                gates = {
                    "exposure": all(e["callable_fraction"] >= .99 for e in estimates.values()),
                    "observed_stopping": all(stops.values()),
                    "two_optimizer_agreement": all(agreement.values()),
                    "target_stability": all(abs(v) < .01 for row in comparisons.values() for v in row["target_relative_change"].values()),
                    "frequency_L1": all(row["frequency_L1"] < .001 for row in comparisons.values()),
                    "segregating_relative_L1": all(row["segregating_relative_L1"] < .05 for row in comparisons.values()),
                    "primary_parameter_recovery": all(abs(v) <= .3 for v in errors["seed1"].values()),
                    "all_estimator_runs_parameter_recovery": all(abs(v) <= .3 for row in errors.values() for v in row.values()),
                    "provenance_bindings": True, "formatter_schema": True,
                }
                primary.extend(abs(v) for v in errors["seed1"].values())
                noise_errors[replicate].extend(abs(v) for v in errors["seed1"].values())
                cohorts[stem] = {"gates": gates, "all_gates_passed": all(gates.values()),
                                 "failed_gates": [key for key, value in gates.items() if not value],
                                 "targets": {label: row["targets"] for label, row in estimates.items()},
                                 "relative_errors": errors, "comparisons": comparisons, "last_epochs": last_epochs}
        assert len(set(seeds)) == 2 and len(primary) == 6
        median = statistics.median(primary)
        passed = all(row["all_gates_passed"] for row in cohorts.values()) and median <= .2
        report = {"verified_at_utc": datetime.now(UTC).isoformat(), "archive_sha256": args.sha256,
                  "freeze_sha256": FREEZE_SHA, "manifest_bound_files_verified": len(manifest["files"]),
                  "cohorts": cohorts, "primary_target_count": len(primary),
                  "primary_median_absolute_relative_error": median,
                  "stratum_medians": {str(k): statistics.median(v) for k, v in noise_errors.items()},
                  "all_registered_gates_passed": passed, "mode": "development", "scientifically_ready_under_frozen_policy": False,
                  "external_binary_bindings": externally_bound,
                  "scope": "Independent raw-spectrum/log/optimizer/target and complete-two-cohort reused-development provenance audit. Raw BAM likelihoods are not rerun; runtime/binary hashes have a separate remote audit.",
                  "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "helper_source_sha256": hashlib.sha256(Path(__file__).with_name("verify_lowcoverage_checkpoint_v2.py").read_bytes()).hexdigest()}
        report["practical_checks"] = {
            name: {"exposure": row["gates"]["exposure"],
                   "optimizer_agreement": row["gates"]["two_optimizer_agreement"],
                   "all_runs_recovery_within_30_percent": row["gates"]["all_estimator_runs_parameter_recovery"],
                   "all_targets_within_5_percent_of_fixed_primary": all(abs(v) <= .05 for c in row["comparisons"].values() for v in c["target_relative_change"].values())}
            for name, row in cohorts.items()}
        report["practically_useful_development_under_separate_policy"] = all(all(row.values()) for row in report["practical_checks"].values())
        report["strict_qualification_passed"] = False
        with args.output.open("x") as handle:
            handle.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"all_registered_gates_passed": passed, "median_absolute_relative_error": median,
                          "failed_gates": {k: v["failed_gates"] for k, v in cohorts.items()}}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    audit(parser.parse_args())
