"""Record final source and completed evidence without altering run artifacts."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

root = Path(__file__).resolve().parent
repo = root.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source_paths = [
    "src/ghist_inference/integration.py",
    "src/ghist_inference/checkpoint.py",
    "src/ghist_inference/freeze.py",
    "src/ghist_inference/validation.py",
    "src/ghist_secondary_contact/spectrum.py",
    "src/ghist_growth_bgs/spectrum.py",
    "src/ghist_split_migration_bgs/spectrum.py",
    "scripts/validate_demography_independently.py",
    "scripts/simulate_growth_bgs_selection.py",
    "scripts/simulate_split_migration_bgs_selection.py",
    "scripts/freeze_secondary_contact.py",
    "scripts/freeze_growth_bgs.py",
    "scripts/freeze_split_migration_bgs.py",
    "scripts/summarize_growth_bgs_fits.py",
    "tests/test_demography_resume.py",
    "tests/test_demography_freeze_cross_review.py",
]
summary = json.loads((root / "summary.json").read_text())
assert all(
    summary[k]["completed"] == summary[k]["planned"]
    for k in ["neutral", "growth_optimizer_status", "neutral_10mb"]
)
assert len(summary["formerly_incomplete"]) == 9
assert set(summary["formerly_incomplete"].values()) == {"completed"}
for name in [
    "comparison_with_previous_completed_runs",
    "growth_optimizer_status_equivalence",
]:
    assert all(
        all(v == 0 for v in model.values())
        for row in summary[name].values()
        for model in row.values()
    )
old_streams, new_streams = set(), []
for directory, target in [("neutral", old_streams), ("neutral_10mb", new_streams)]:
    for path in (root / directory).glob("task*.json"):
        row = json.loads(path.read_text())
        if "configuration" not in row:
            continue
        if isinstance(target, set):
            target.update(row["seeds"].values())
        else:
            target.extend(row["seeds"].values())
assert len(new_streams) == len(set(new_streams)) == 112
assert not old_streams.intersection(new_streams)
checks = json.loads((root / "likelihood_refinement.json").read_text())["checks"]
checks = [record for models in checks.values() for record in models.values()]
assert len(checks) == 53
linux_base = "/home/kkor/ghist-2026-review-20261120"
verification = {
    "date": "2026-11-20",
    "base_commit": "b847c4c",
    "final_data_runs_performed": False,
    "submissions": 0,
    "public_answers_used": False,
    "scope": "Numerical repairs, complete diagnostic recovery and explicit scientific failures; not a readiness certification",
    "completed_datasets": {
        "neutral_2mb": 22,
        "neutral_10mb": 28,
        "same_spectrum_growth_diagnostic_replays": 8,
        "previous_unscaled_split_tree_resumed": 1,
    },
    "censored_or_failed_datasets": 0,
    "selected_nonconverged_fits": {
        "neutral_2mb": 4,
        "neutral_10mb": 7,
        "termination": "L-BFGS-B warnflag=2, ABNORMAL; retained and rejected by freeze guard",
    },
    "streams": {
        "new_panel_streams": 112,
        "unique_new_panel_streams": 112,
        "overlap_with_specific_2mb_panel": 0,
        "legacy_scheme_does_reuse_fit0_fit1": True,
    },
    "tests": {
        "command": ".venv/bin/python -m pytest -q tests/test_demography_freeze_cross_review.py tests/test_demography_resume.py tests/test_secondary_contact_scientific.py tests/test_secondary_contact.py tests/test_bgs_physical_units.py tests/test_expected_spectrum_domain.py tests/test_reliability.py tests/test_growth_bgs.py tests/test_split_migration_bgs.py",
        "passed": 128,
        "warnings": 1,
        "warning": "PyTorch nested-tensor performance warning",
        "critical_ruff": "passed",
    },
    "cross_review": {
        "added_regression_cases": 18,
        "findings_fixed": [
            "Duplicate normalized growth projections could hide an unconverged fit while contributing to the reported median; duplicate inputs are now rejected before aggregation.",
            "Incomplete converged projection sets could freeze; Task3 now requires exactly 84,60,40 and Task4 exactly 50x40,70x50,94x68.",
            "Task2 real freeze main could export a nonfinite final prediction; final, raw, refined and nominal risk values now require finite physical domains, and strict atomic JSON rejects nonfinite residual fields.",
            "Unknown fit checkpoint statuses could silently call the fitter and overwrite the record; they are now rejected before callback or write.",
        ],
        "valid_boundaries_preserved": [
            "Task2 contact time zero or equal to split time",
            "Task2 migration rate zero",
            "An unchanged zero-valued refinement reference",
        ],
        "numerical_model_equations_changed": False,
        "historical_run_artifacts_changed": False,
        "no_additional_scientific_readiness_claim": True,
    },
    "runtime": {
        "local": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": {
                k: importlib.metadata.version(k)
                for k in ["moments-popgen", "numpy", "scipy", "msprime", "tskit"]
            },
        },
        "linux": json.loads((root / "linux_runtime.json").read_text()),
        "betty": json.loads((root / "betty_runtime.json").read_text()),
    },
    "commands": {
        "neutral_2mb": f"cd {linux_base}; LD_LIBRARY_PATH=$PWD/.native/usr/lib/aarch64-linux-gnu PYTHONPATH=$PWD/src:$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-reliability/bin/python -m scripts.validate_demography_independently --out artifacts/readiness/demography/neutral --workers 7 --replicates 2 --restarts 12 --seed-scheme legacy_v1",
        "growth_diagnostic_replay": f"cd {linux_base}; LD_LIBRARY_PATH=$PWD/.native/usr/lib/aarch64-linux-gnu PYTHONPATH=$PWD/src:$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-reliability/bin/python -m scripts.validate_demography_independently --out artifacts/readiness/demography/growth_optimizer_status --workers 1 --replicates 2 --restarts 12 --seed-scheme legacy_v1 --tasks 3",
        "growth_replay_preparation": "Copy the eight neutral/task3*.sfs.npz to growth_optimizer_status/ before running; no fit checkpoints copied.",
        "neutral_10mb": f"cd {linux_base}/snapshots/readiness-demography-final; LD_LIBRARY_PATH={linux_base}/.native/usr/lib/aarch64-linux-gnu PYTHONPATH=$PWD/src:$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 {linux_base}/.venv-reliability/bin/python -m scripts.validate_demography_independently --out {linux_base}/artifacts/readiness/demography/neutral_10mb --workers 4 --replicates 4 --restarts 12 --sequence-length 10000000 --seed-scheme counter_v2 --tasks 3 4",
        "unscaled_split": "sbatch experiments/readiness_2026_11_20/demography/resume_split_q1.sbatch (Betty paths are encoded in the recorded script)",
        "native_reproducer": ".venv/bin/python experiments/readiness_2026_11_20/demography/reproduce_native_nonprogress.py",
        "summary": ".venv/bin/python experiments/readiness_2026_11_20/demography/summarize.py",
        "timestep_check": ".venv/bin/python experiments/readiness_2026_11_20/demography/check_likelihood_refinement.py",
    },
    "likelihood_refinement": {
        "fitted_models": 53,
        "maximum_absolute_log_likelihood_change": max(
            abs(r["log_likelihood_change"]) for r in checks
        ),
        "maximum_original_score_recomputation_difference": max(
            abs(r["recorded_score_difference"]) for r in checks
        ),
        "maximum_relative_change_positive_observed_bin": max(
            r["maximum_relative_change_in_positive_observed_bins"] for r in checks
        ),
        "all_observed_bins_retained": all(
            r["all_observed_bins_included"] for r in checks
        ),
    },
    "source_sha256_final": {p: digest(repo / p) for p in source_paths},
    "source_snapshot_caveats": [
        "The 2Mb panel began before final input-control validation and exact-bound initializer relaxation; generic starts and all selected inputs were unchanged.",
        "The original 2Mb Task3 fit records predate optimizer metadata; separate same-spectrum diagnostic replays match all estimates and likelihoods exactly.",
        "The immutable 10Mb snapshot includes final numerical source. Subsequent changes add export guards, preserve growth summary diagnostics, retry failed benchmark checkpoints and correct descriptive scope; they do not change successful model fitting.",
        "The 10Mb raw run scope string incorrectly says 2Mb; the generating length field and checkpoint metadata are 10000000. Raw records are preserved.",
        "profile_censored.py and profile_censored_bounded.py target baseline b847c4c and its five-item run_one tuple; reproduce_native_nonprogress.py is the current standalone reproducer.",
        "Failed initial unbounded trace attempts provide no completed-run evidence; only the bounded independent reproducer establishes the native nonprogress finding.",
        "Final source hashes are distinct from runtime snapshot hashes and are not retroactively attributed to earlier processes.",
        "Independent cross-review subsequently added duplicate/exact projection checks, finite physical-domain Task2 export checks and unknown-status checkpoint rejection. These guards do not change successful model fitting and are not attributed to earlier runtime snapshots.",
    ],
    "artifacts": {
        str(p.relative_to(root)): {"sha256": digest(p), "bytes": p.stat().st_size}
        for p in sorted(root.rglob("*"))
        if p.is_file()
        and p.name != "verification.json"
        and "__pycache__" not in p.parts
    },
}
(root / "verification.json").write_text(
    json.dumps(verification, indent=2, allow_nan=False) + "\n"
)
print(
    "Verified 22+28 completed neutral datasets, 8 identical-spectrum replays, 1 resumed Q1 BGS fit, 53 timestep checks."
)
