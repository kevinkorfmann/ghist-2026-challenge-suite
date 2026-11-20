"""Fixed full-cohort prospective pragmatic assessment; retain strict diagnostics."""

from pathlib import Path
import argparse
import itertools
import json
from calibrate_phase import metrics
from structural_confidence import reliability
from assess_structural_qualification import numerical_gates
from ghist_relatedness.labels import HUMAN_LABELS
from scripts.score_relatedness_published_metric import evaluate_published_metric
from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--accounting", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    verify_inputs(json.loads((source / "source_manifest.json").read_text()))
    for name in ["detected-panel-source", "structural-qualification-source-v3"]:
        verify_inputs(
            json.loads((args.root.parent / name / "source_manifest.json").read_text())
        )
    accounting = json.loads(args.accounting.read_text())
    labels = list(HUMAN_LABELS.values())
    variants = {name: [] for name in ["structural", "pair_only", "baseline"]}
    bound = {str(args.accounting): sha256(args.accounting)}
    failures = []
    resources = {}
    methods = {}
    for index in [8, 9]:
        folder = args.root / f"cohort_{index:03d}"
        required = [
            folder / name
            for name in [
                "result.json",
                "truth.json",
                "identifiers.json",
                "structural_predictions.json",
                "pair_only_predictions.json",
                "baseline_predictions.json",
                "structural_tree.json",
                "missingness_contract/result.json",
            ]
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            failures.append({"index": index, "missing": missing})
            continue
        bound.update({str(path): sha256(path) for path in required})
        record = json.loads((folder / "result.json").read_text())
        verify_inputs(record["output_sha256"])
        roster = sorted(json.loads((folder / "identifiers.json").read_text()).values())
        expected = set(itertools.combinations(roster, 2))
        truth = json.loads((folder / "truth.json").read_text())
        if (
            len(roster) != 229
            or len(expected) != 26106
            or len(truth) != 26106
            or {tuple(sorted([row["indID1"], row["indID2"]])) for row in truth}
            != expected
        ):
            raise ValueError("Incomplete fresh229roster or truth")
        for name in variants:
            predictions = json.loads((folder / f"{name}_predictions.json").read_text())
            if (
                len(predictions) != 26106
                or {
                    tuple(sorted([row["indID1"], row["indID2"]])) for row in predictions
                }
                != expected
            ):
                raise ValueError("Incomplete fresh prediction set")
            variants[name].append(
                {
                    "cohort_id": f"nested_{index:03d}",
                    "truth": truth,
                    "predictions": predictions,
                }
            )
        audit = accounting.get(str(index), {})
        resources[str(index)] = {
            "completed": audit.get("state") == "COMPLETED",
            "matching_job": audit.get("job_id") == record["slurm_job_id"],
            "wall_time": 0 < audit.get("elapsed_seconds", 0) <= 5400,
            "memory": 0 < audit.get("max_rss_kib", 0) <= 44 * 1024**2,
        }
        contract = json.loads((folder / "missingness_contract/result.json").read_text())
        if (
            contract["state"] != "MISSINGNESS_REJECTION_VERIFIED"
            or contract["biological_robustness_qualified"]
        ):
            failures.append(
                {"index": index, "reason": "Missingness rejection contract failed"}
            )
        methods[str(index)] = json.loads((folder / "structural_tree.json").read_text())[
            "pragmatic_policy"
        ]
    args.output.mkdir(parents=True, exist_ok=False)
    for name, cohorts in variants.items():
        save_json(
            args.output / f"{name}_scorer_input.json",
            {"labels": labels, "cohorts": cohorts},
        )
    if failures or len(variants["structural"]) != 2:
        save_json(
            args.output / "result.json",
            {
                "state": "INCOMPLETE_PROSPECTIVE_PRAGMATIC_SCALING_FAILED",
                "passed": False,
                "failures": failures,
                "resources": resources,
                "missing_cohorts_excluded": False,
                "input_sha256": bound,
                "strict_v3_status": "FAILED",
            },
        )
        return
    scores = {}
    published = {}
    strict = {}
    pragmatic = {}
    reliability_reports = {}
    for label, selection in [
        ("cohort_008", [0]),
        ("cohort_009", [1]),
        ("pooled", [0, 1]),
    ]:
        selected = {
            name: [cohorts[i] for i in selection] for name, cohorts in variants.items()
        }
        scores[label] = {
            name: metrics(cohorts, labels) for name, cohorts in selected.items()
        }
        published[label] = {
            name: evaluate_published_metric(cohorts, labels)
            for name, cohorts in selected.items()
        }
        reliability_reports[label] = reliability(selected["structural"], labels)
        current = scores[label]["structural"]
        base = scores[label]["baseline"]
        strict[label] = numerical_gates(current, base) | {
            "ece": reliability_reports[label]["all"]["ece_10bins"] <= 0.1
        }
        pragmatic[label] = {
            "accuracy": current["exact_accuracy"] >= 0.70,
            "macro_recall": current["macro_recall"] >= 0.60,
            "weighted_jaccard_exceeds_baseline": published[label]["structural"][
                "pooled"
            ]["weighted_jaccard"]
            > published[label]["baseline"]["pooled"]["weighted_jaccard"],
            "complete_all_pairs": published[label]["structural"][
                "pair_coverage_passed"
            ],
        }
    pragmatic["resources"] = {
        "all": all(all(values.values()) for values in resources.values())
    }
    verify_inputs(bound)
    save_json(
        args.output / "result.json",
        {
            "state": "ASSESSED_NEW_POLICY_PROSPECTIVE_SCALING",
            "passed": all(all(values.values()) for values in pragmatic.values()),
            "pragmatic_gates": pragmatic,
            "strict_numeric_diagnostics": strict,
            "strict_v3_status": "FAILED_UNCHANGED",
            "scores": scores,
            "published_metric": published,
            "reliability": reliability_reports,
            "resources": resources,
            "selected_method": methods,
            "input_sha256": bound,
            "source_manifest_sha256": sha256(source / "source_manifest.json"),
            "scope": "Two fresh229person complete-call cohorts under fixed nested generator; final input263people and structured missingness remain extrapolations; ECE not a pragmatic pass requirement; no submission performed",
        },
    )


if __name__ == "__main__":
    main()
