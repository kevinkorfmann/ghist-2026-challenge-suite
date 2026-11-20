"""Run one prespecified biological case through the production VCF predictor.

The separate prediction process receives observations and accessibility only.
Truth is consumed for scoring only after candidate files have been written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import warnings
from pathlib import Path

import msprime
import stdpopsim
import tskit

from ghist_single_sweep.export import write_biallelic_vcf

from ghist_single_sweep.peaks import highest_peak_position, read_raisd
from ghist_single_sweep.scan import read_bed_intervals
from ghist_submission.validate import validate_submission


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def model_for(name):
    if name == "Zigzag_1S14":
        return stdpopsim.get_species("HomSap").get_demographic_model(name), "generic"
    if name == "constant_10000":
        demography = msprime.Demography()
        demography.add_population(name="generic", initial_size=10_000)
        return stdpopsim.DemographicModel(
            id="ReadinessConstant",
            description="Constant diploid N=10000",
            long_description="Prespecified stress history; not an asserted challenge history.",
            generation_time=25,
            model=demography,
            mutation_rate=1.44e-8,
        ), "generic"
    raise ValueError(f"Unknown prespecified history: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--case", type=int, required=True)
    parser.add_argument("--accessible", type=Path, required=True)
    parser.add_argument("--slim", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--reuse-case-dir", type=Path)
    args = parser.parse_args()
    design = json.loads(args.design.read_text())
    case = design["cases"][args.case]
    folder = args.out_dir / f"case_{args.case:03d}"
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    status = {
        "case": args.case,
        "seed": case["seed"],
        "state": "RUNNING",
        "design_sha256": sha256(args.design),
    }
    (folder / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    try:
        if sha256(args.accessible) != design["accessibility_sha256"]:
            raise ValueError("Accessibility differs from the frozen design")
        for source, digest in design["source_sha256"].items():
            if sha256(Path(source)) != digest:
                raise ValueError(f"Source differs from frozen design: {source}")
        if args.reuse_case_dir:
            previous = json.loads((args.reuse_case_dir / "truth.json").read_text())
            if any(previous[key] != value for key, value in case.items()):
                raise ValueError("Reused simulation differs from the frozen case")
            if (
                previous["slim_scaling_factor"] != 1
                or previous["slim_burn_in"] != design["generation"]["slim_burn_in"]
                or previous["mutation_rate"] != 1.44e-8
            ):
                raise ValueError(
                    "Reused simulation uses different physical/scaling settings"
                )
            ts = tskit.load(args.reuse_case_dir / "simulation.trees")
            status["reused_tree_sha256"] = sha256(
                args.reuse_case_dir / "simulation.trees"
            )
            status["reused_resource_pilot_biology"] = str(args.reuse_case_dir)
        else:
            model, population = model_for(case["history"])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                contig = stdpopsim.get_species("HomSap").get_contig(
                    "7",
                    genetic_map="HapMapII_GRCh38",
                    mutation_rate=1.44e-8,
                    # stdpopsim 0.3.0's BED reader squeezes a single row to a scalar.
                    # Our validated reader preserves the required (n, 2) shape.
                    inclusion_mask=read_bed_intervals(args.accessible, "7"),
                )
            if any("Genetic map not found" in str(item.message) for item in caught):
                raise ValueError("Missing genetic map: flat-map fallback forbidden")
            contig.add_single_site(id="target", coordinate=case["position"])
            sweep = {
                "single_site_id": "target",
                "population": population,
                "mutation_generation_ago": case["origin_generations"],
                "selection_coeff": case["selection_coefficient"],
                "min_freq_at_end": 0.95,
            }
            if "selection_start_generations" in case:
                sweep["start_generation_ago"] = case["selection_start_generations"]
                sweep["min_freq_at_start"] = case["minimum_start_frequency"]
            ts = stdpopsim.get_engine("slim").simulate(
                model,
                contig,
                samples={population: 25},
                seed=case["seed"],
                extended_events=stdpopsim.selective_sweep(**sweep),
                slim_path=str(args.slim),
                slim_scaling_factor=1,
                slim_burn_in=design["generation"]["slim_burn_in"],
            )
        if ts.num_samples != 50:
            raise ValueError(
                f"Expected 25 diploids / 50 haplotypes, got {ts.num_samples}"
            )
        ts.dump(folder / "simulation.trees")
        truth = {
            **case,
            "diploid_samples": 25,
            "haplotypes": ts.num_samples,
            "sequence_length": ts.sequence_length,
            "stdpopsim_version": stdpopsim.__version__,
            "mutation_rate": 1.44e-8,
            "background_selection": False,
            "slim_scaling_factor": 1,
            "slim_burn_in": design["generation"]["slim_burn_in"],
            "conditioning": "SLiM rejection conditions selected allele frequency >=.95 at sampling; establishment probability not estimated",
        }
        (folder / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
        vcf = folder / "simulation.vcf"
        export = write_biallelic_vcf(ts, vcf, "7")
        if export["contig_length"] != design["generation"]["length"]:
            raise ValueError("Exported contig length differs from target exposure")
        (folder / "export.json").write_text(json.dumps(export, indent=2) + "\n")
        del ts
        status["simulation_seconds"] = time.monotonic() - started
        with (folder / "prediction.log").open("w") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "scripts.predict_single_sweep",
                    "--vcf",
                    str(vcf),
                    "--accessible",
                    str(args.accessible),
                    "--dataset",
                    "testing",
                    "--bed",
                    str(folder / "candidate.bed"),
                    "--out",
                    str(folder / "prediction.json"),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        status["prediction_seconds"] = (
            time.monotonic() - started - status["simulation_seconds"]
        )
        # Mechanistically distinct whole-chromosome reference, without truth.
        with (folder / "reference.log").open("w") as log:
            subprocess.run(
                [
                    str(args.raisd.resolve()),
                    "-n",
                    "reference",
                    "-I",
                    str(vcf.resolve()),
                    "-w",
                    "50",
                    "-R",
                    "-f",
                ],
                cwd=folder,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        report = read_raisd(folder / "RAiSD_Report.reference")
        reference_info = (folder / "RAiSD_Info.reference").read_text()
        reference_sites = [
            int(value) for value in re.findall(r"\|\s+snps\s+(\d+)", reference_info)
        ]
        if reference_sites != [export["accepted_biallelic_polymorphic_sites"]]:
            raise ValueError(
                "RAiSD and production predictor did not use identical SNP sets"
            )
        status["reference_usable_sites"] = reference_sites[0]
        reference_position = highest_peak_position(report, export["contig_length"])
        prediction = json.loads((folder / "prediction.json").read_text())
        if (
            prediction["contig_length"] != export["contig_length"]
            or prediction["haplotypes"] != 50
            or prediction["accepted_biallelic_sites"]
            != export["accepted_biallelic_polymorphic_sites"]
        ):
            raise ValueError(
                "Production reader changed exported site/exposure/sample contract"
            )
        validation = validate_submission(folder / "candidate.bed", 5, "testing")
        lo, hi = prediction["interval"]
        metrics = {
            "case": args.case,
            "seed": case["seed"],
            "history": case["history"],
            "selection_coefficient": case["selection_coefficient"],
            "origin_generations": case["origin_generations"],
            "method": prediction["method"],
            "consensus_eligible": prediction["consensus_eligible"],
            "absolute_error": abs(prediction["center"] - case["position"]),
            "covered": lo <= case["position"] < hi,
            "interval_width": hi - lo,
            "format_validation": validation,
            "reference_position": reference_position,
            "reference_absolute_error": abs(reference_position - case["position"]),
            "reference_report_sha256": sha256(folder / "RAiSD_Report.reference"),
            "reference_info_sha256": sha256(folder / "RAiSD_Info.reference"),
            "prediction_sha256": sha256(folder / "prediction.json"),
            "bed_sha256": sha256(folder / "candidate.bed"),
        }
        (folder / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
        status["state"] = "COMPLETED"
    except Exception as error:
        status.update(state="FAILED", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        status["elapsed_seconds"] = time.monotonic() - started
        (folder / "status.json").write_text(json.dumps(status, indent=2) + "\n")
        print(json.dumps(status))


if __name__ == "__main__":
    main()
