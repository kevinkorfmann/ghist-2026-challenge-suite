"""Recompute descriptive numerical comparisons from all 48 planned fits."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    directory = ROOT / "results"
    protocol = json.loads((directory / "protocol.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    if digest(directory / "protocol.json") != completion["protocol_sha256"]:
        raise ValueError("Protocol hash mismatch")
    for name, expected in completion["output_sha256"].items():
        if digest(directory / name) != expected:
            raise ValueError(f"Raw output hash mismatch: {name}")
    if digest(ROOT / "diagnose.py") != protocol["source_sha256"]["diagnostic"]["sha256"]:
        raise ValueError("Diagnostic source hash mismatch")
    for filename, expected in protocol["input_sha256"].items():
        name = Path(filename).name
        path = (ROOT.parent / "results/protocol.json") if name == "panel_protocol.json" else ROOT.parent / "results/case_08" / name
        if digest(path) != expected:
            raise ValueError(f"Original local input changed: {name}")
    rows = json.loads((directory / "fits.json").read_text())
    original = json.loads((directory / "original_selected_diagnostics.json").read_text())
    expected = {(method, variant["name"], index)
                for method in protocol["methods"] for variant in protocol["variants"]
                for index in range(protocol["restarts_per_variant"])}
    identities = [(row["method"], row["variant"], row["restart"]) for row in rows]
    if len(identities) != len(set(identities)) or set(identities) != expected:
        raise ValueError("Missing, duplicate or unexpected fit records")
    groups = []
    for method in protocol["methods"]:
        for variant in protocol["variants"]:
            group = [row for row in rows if row["method"] == method and row["variant"] == variant["name"]]
            finite = [row for row in group if "objective" in row]
            best = min(finite, key=lambda row: row["objective"])
            groups.append({
                "method": method, "variant": variant["name"],
                "planned_restarts": protocol["restarts_per_variant"], "finite_restarts": len(finite),
                "reported_converged": sum(row.get("optimizer", {}).get("converged", False) for row in group),
                "best_restart": best["restart"], "best_parameters": best["parameters"],
                "best_objective": best["objective"], "best_reported_converged": best["optimizer"]["converged"],
                "best_termination": best["optimizer"]["message"],
                "best_centered_gradient_inf_norm_h1e-5": max(abs(item["central_slope"]) for item in best["local_differences"]["1e-05"]),
                "worst_finite_objective": max(row["objective"] for row in finite),
                "restarts_within_0_1_of_best": sum(row["objective"] - best["objective"] <= 0.1 for row in finite),
            })
    summary = {
        "scope": "Numerical diagnostic selected for abnormal termination; no independent biological recovery or promotion claim",
        "job_id": "10851187", "elapsed_seconds": 107, "allocated_cpus": 1, "max_rss_kb": 236788,
        "planned_fits": len(expected), "recorded_fits": len(rows), "censored_fits": 0,
        "all_local_original_input_hashes_verified": True,
        "raw_output_hashes_verified": True,
        "groups": groups,
        "original_selected": {method: {
            "recomputed_minus_stored_objective": row["recomputed_objective"] - row["stored_objective"],
            "coarse_forward_slopes": [item["forward_slope"] for item in row["local_differences"]["0.001"]],
            "fine_centered_slopes": [item["central_slope"] for item in row["local_differences"]["1e-05"]],
        } for method, row in original.items()},
        "truth_used": False, "original_panel_modified": False, "production_source_modified": False,
    }
    with tempfile.NamedTemporaryFile(mode="w", dir=ROOT, delete=False) as handle:
        handle.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    os.replace(handle.name, ROOT / "summary.json")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
