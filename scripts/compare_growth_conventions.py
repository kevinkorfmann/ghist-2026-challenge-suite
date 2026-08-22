from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_single(path: Path) -> tuple[int, dict]:
    payload = json.loads(path.read_text())
    if len(payload["results"]) != 1:
        raise ValueError(f"{path} must contain exactly one fitted stratum")
    return int(payload["projection"]), next(iter(payload["results"].values()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exponential", type=Path, action="append", required=True)
    parser.add_argument("--sudden", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    exponential = dict(load_single(path) for path in args.exponential)
    sudden = dict(load_single(path) for path in args.sudden)
    if set(exponential) != set(sudden):
        raise ValueError("growth conventions must have identical projections")

    comparisons = {}
    for projection in sorted(exponential):
        exp_fit = exponential[projection]
        sudden_fit = sudden[projection]
        exp_ll = exp_fit["multinomial_log_likelihood"]
        sudden_ll = sudden_fit["multinomial_log_likelihood"]
        comparisons[str(projection)] = {
            "exponential": {
                "generations": exp_fit["generations"],
                "post_growth_scaling": exp_fit["post_growth_scaling"],
                "multinomial_log_likelihood": exp_ll,
            },
            "sudden": {
                "generations": sudden_fit["generations"],
                "post_growth_scaling": sudden_fit["post_growth_scaling"],
                "multinomial_log_likelihood": sudden_ll,
            },
            "log_likelihood_exponential_minus_sudden": exp_ll - sudden_ll,
            "preferred": "exponential" if exp_ll > sudden_ll else "sudden",
        }
    payload = {
        "stratum": "noncoding_1kb_100kb",
        "equal_parameter_count": True,
        "comparisons": comparisons,
        "unanimous_preference": len(
            {row["preferred"] for row in comparisons.values()}
        )
        == 1,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
