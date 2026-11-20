"""Check best-fit scores and rare observed bins at half the fitted timestep."""

import json
from pathlib import Path
import moments
import numpy as np
from ghist_inference.validation import validate_expected_spectrum
from ghist_secondary_contact.spectrum import (
    secondary_contact_spectrum,
    generator_matched_secondary_contact_spectrum,
)
from ghist_split_migration_bgs.spectrum import isolation_with_migration_spectrum

root = Path(__file__).resolve().parent
checks = {}


def check(spectrum, fit, function, parameters):
    data = moments.Spectrum(spectrum).project(fit["projection"])

    def evaluate(step):
        expected = function(parameters, data.sample_sizes, dt_fac=step, adapt_dt=False)
        try:
            return validate_expected_spectrum(expected, data)
        except ValueError:
            return validate_expected_spectrum(
                function(parameters, data.sample_sizes, dt_fac=step, adapt_dt=True),
                data,
            )

    coarse = evaluate(fit["integration_dt_fac"])
    fine = evaluate(fit["integration_dt_fac"] / 2)
    active = ~np.ma.getmaskarray(data)
    positive = active & (np.asarray(data) > 0)
    coarse_values, fine_values = np.asarray(coarse), np.asarray(fine)
    coarse_ll = float(moments.Inference.ll_multinom(coarse, data))
    fine_ll = float(moments.Inference.ll_multinom(fine, data))
    return {
        "recorded_log_likelihood": fit["multinomial_log_likelihood"],
        "recomputed_log_likelihood": coarse_ll,
        "halved_timestep_log_likelihood": fine_ll,
        "log_likelihood_change": fine_ll - coarse_ll,
        "recorded_score_difference": coarse_ll - fit["multinomial_log_likelihood"],
        "relative_L1_spectrum_change": float(
            np.abs(fine_values[active] - coarse_values[active]).sum()
            / fine_values[active].sum()
        ),
        "maximum_relative_change_in_positive_observed_bins": float(
            np.max(
                np.abs(fine_values[positive] - coarse_values[positive])
                / fine_values[positive]
            )
        ),
        "all_observed_bins_included": True,
    }


for path in sorted(
    [
        *(root / "neutral").glob("task*.json"),
        *(root / "neutral_10mb").glob("task*.json"),
    ]
):
    row = json.loads(path.read_text())
    if "configuration" not in row or row["configuration"]["task"] == 3:
        continue
    spectrum = np.load(path.with_suffix(".sfs.npz"))["spectrum"]
    checks[str(path.relative_to(root))] = {}
    for name, fit in row["fits"].items():
        p = fit["scaled_parameters"]
        if row["configuration"]["task"] == 2:
            parameters = [p["population1_size_ratio"], p["population2_size_ratio"]]
            if name == "simple":
                function = secondary_contact_spectrum
                parameters += [
                    p["isolation_duration_2N"],
                    p["contact_duration_2N"],
                    p["migration_2N"],
                ]
            else:
                function = generator_matched_secondary_contact_spectrum
                parameters += [
                    p["split_duration_2N"],
                    p["contact_fraction"],
                    p["migration_2N"],
                ]
        else:
            function = isolation_with_migration_spectrum
            parameters = [
                p[key]
                for key in (
                    "AFR_size_ratio",
                    "EUR_size_ratio",
                    "split_duration_2N",
                    "migration_2N",
                )
            ]
        checks[str(path.relative_to(root))][name] = check(
            spectrum, fit, function, parameters
        )
path = root / "bgs/split_q1.json"
if path.exists():
    row = json.loads(path.read_text())
    fit = row["fit"]
    p = fit["scaled_parameters"]
    parameters = [
        p[key]
        for key in (
            "AFR_size_ratio",
            "EUR_size_ratio",
            "split_duration_2N",
            "migration_2N",
        )
    ]
    spectrum = np.load(path.with_suffix(".sfs.npz"))["spectrum"]
    checks["split_q1.json"] = {
        "fit": check(spectrum, fit, isolation_with_migration_spectrum, parameters)
    }
(root / "likelihood_refinement.json").write_text(
    json.dumps(
        {
            "scope": "Fixed fitted parameters, full observed likelihood domain; timestep sensitivity only, not accuracy or global optimizer proof",
            "checks": checks,
        },
        indent=2,
    )
    + "\n"
)
print(f"Checked {sum(map(len, checks.values()))} fitted models")
