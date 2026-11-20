"""Verify source-only follow-up and preservation of prior experiment evidence."""

from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import random
import subprocess
import sys
import xml.etree.ElementTree as ET

DIRECTORY = Path(__file__).resolve().parent
REPO = DIRECTORY.parents[2]
sys.path[:0] = [str(REPO / "src"), str(REPO)]

from scripts.consensus_task11_dog_pedigrees import METHOD_VERSION  # noqa: E402
from scripts.reconstruct_bonsai_pedigree import coverage_intervals  # noqa: E402


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    preserved = json.loads((DIRECTORY / "preserved_artifacts_before.json").read_text())
    assert len(preserved) == 171
    assert all(digest(REPO / path) == expected for path, expected in preserved.items())
    before = json.loads((DIRECTORY / "source_before.json").read_text())
    after = {path: digest(REPO / path) for path in before}
    junit = ET.parse(DIRECTORY / "pytest.xml").getroot()
    suites = list(junit.iter("testsuite"))
    counts = {
        key: sum(int(suite.get(key, "0")) for suite in suites)
        for key in ("tests", "errors", "failures", "skipped")
    }
    assert counts == {"tests": 170, "errors": 0, "failures": 0, "skipped": 0}

    # Pointwise active-edge oracle; independent of the converter's event sweep.
    rng = random.Random(2026112009)
    for _ in range(5000):
        segments = []
        for _ in range(rng.randrange(0, 25)):
            start, end = sorted(rng.sample(range(1, 31), 2))
            segments.append((start, end, rng.choice((1, 2)), rng.choice((1, 2))))
        result = coverage_intervals(segments)
        for position in range(1, 30):
            active = {
                (left, right)
                for start, end, left, right in segments
                if start <= position < end
            }
            full = any(
                left[0] != right[0] and left[1] != right[1]
                for left, right in combinations(active, 2)
            )
            observed = [
                state for start, end, state in result if start <= position < end
            ]
            assert observed == ([full] if active else [])

    lint_paths = [
        "scripts/consensus_task11_dog_pedigrees.py",
        "scripts/prepare_hapibd_vcf.py",
        "tests/test_task11_consensus.py",
        "tests/test_hapibd_preparation.py",
        str(Path(__file__).resolve().relative_to(REPO)),
    ]
    for command in (
        [str(REPO / ".venv/bin/ruff"), "check", *lint_paths],
        ["git", "diff", "--check"],
    ):
        subprocess.run(command, cwd=REPO, check=True, capture_output=True, text=True)
    record = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Source contract cross-review; no full-pipeline rerun or scientific promotion",
        "method_version": METHOD_VERSION,
        "findings_fixed": [
            "Exact ties are independent of probability mapping order",
            "Only explicit (+inf,+inf,None) unrelated sentinel bypasses finite topology validation",
            "Retained heterozygotes cannot mix or cross independent phase sets per sample/chromosome",
        ],
        "tests": {
            "command": ".venv/bin/python -m pytest -q tests/test_bonsai_conversion.py tests/test_hapibd_preparation.py tests/test_task11_consensus.py tests/test_task11_freeze_guard.py tests/test_relatedness_missingness.py tests/test_task10_complete_pipeline.py tests/test_task11_complete_pipeline.py --junitxml=experiments/complete_pipeline_2026_11_20/cross_review/pytest.xml",
            **counts,
            "junit_sha256": digest(DIRECTORY / "pytest.xml"),
        },
        "independent_geometry_check": {
            "seed": 2026112009,
            "random_interval_sets": 5000,
            "point_checks_per_set": 29,
            "passed": True,
            "scope": "Synthetic software geometry, not independent biological validation",
        },
        "ruff": "passed",
        "git_diff_check": "passed",
        "preserved_experiment_files": len(preserved),
        "preserved_manifest_sha256": digest(
            DIRECTORY / "preserved_artifacts_before.json"
        ),
        "all_prior_experiment_hashes_unchanged": True,
        "source_sha256_before": before,
        "source_sha256_after": after,
        "verifier_sha256": digest(Path(__file__)),
        "new_scientific_runs": 0,
        "historical_results_rewritten": False,
        "final_data_runs_performed": False,
        "submissions": 0,
        "public_answers_used": False,
        "runtime_attribution": "Earlier followup_verification.json and runtime snapshots retain their original source hashes; they do not certify these later source changes.",
    }
    (DIRECTORY / "verification.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n"
    )
    print("Verified 170 tests, 5000 interval sets and 171 unchanged experiment files.")


if __name__ == "__main__":
    main()
