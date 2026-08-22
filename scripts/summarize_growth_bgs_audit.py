from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stratum_summary(source: dict) -> dict:
    sites = int(source["biallelic_sites"])
    folded = source["folded_sfs"]
    return {
        "segregating_biallelic_sites": sites,
        "oriented_sites": int(source["oriented_sites"]),
        "singleton_fraction": int(folded.get("1", 0)) / max(sites, 1),
        "doubleton_fraction": int(folded.get("2", 0)) / max(sites, 1),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic-audit", type=Path, required=True)
    parser.add_argument("--resource-audit", type=Path, required=True)
    parser.add_argument("--generic-job", required=True)
    parser.add_argument("--resource-job", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    generic = json.loads(args.generic_audit.read_text())
    resource = json.loads(args.resource_audit.read_text())
    windows = resource["windows_100kb"]
    site_counts = np.asarray([row["biallelic_sites"] for row in windows])
    nonempty = site_counts > 0
    diversity = np.asarray([row["diversity_per_bp"] for row in windows])
    coding = np.asarray([row["coding_fraction"] for row in windows])
    recombination = np.asarray(
        [row["mean_variant_site_recombination_cm_per_mb"] for row in windows]
    )
    by_contig = {}
    for contig in resource["contig_lengths"]:
        mask = np.asarray([row["contig"] == contig for row in windows]) & nonempty
        by_contig[contig] = {
            "nonempty_windows": int(mask.sum()),
            "mean_diversity_per_bp": float(diversity[mask].mean()),
            "median_diversity_per_bp": float(np.median(diversity[mask])),
        }
    result = {
        "jobs": {
            "generic_audit": args.generic_job,
            "resource_audit": args.resource_job,
        },
        "source_audits": {
            "generic": {
                "path": str(args.generic_audit),
                "sha256": sha256(args.generic_audit),
            },
            "resource": {
                "path": str(args.resource_audit),
                "sha256": sha256(args.resource_audit),
            },
        },
        "files": resource["files"],
        "samples": resource["samples"],
        "expected_haplotypes": resource["expected_haplotypes"],
        "sequence_length": resource["sequence_length"],
        "contig_lengths": resource["contig_lengths"],
        "records": resource["records"],
        "biallelic_sites": resource["biallelic_sites"],
        "segregating_biallelic_sites": resource["segregating_biallelic_sites"],
        "nonsegregating_biallelic_sites": resource[
            "nonsegregating_biallelic_sites"
        ],
        "multiallelic_sites": resource["multiallelic_sites"],
        "ancestral_orientation": resource["ancestral_orientation"],
        "phased_call_fraction": generic["phased_call_fraction"],
        "sample_missingness_quantiles": generic["sample_missingness_quantiles"],
        "vcf_info_fields": resource["vcf_info_fields"],
        "vcf_format_fields": resource["vcf_format_fields"],
        "nonsynonymous_annotation_present": resource[
            "nonsynonymous_annotation_present"
        ],
        "coding_bed": resource["coding_bed"],
        "coding_union_bp": sum(
            record["union_bp"] for record in resource["coding_bed"].values()
        ),
        "coding_fraction": sum(
            record["union_bp"] for record in resource["coding_bed"].values()
        )
        / resource["sequence_length"],
        "recombination_map": resource["recombination_map"],
        "positive_map_rate_quartiles_cm_per_mb": resource[
            "positive_map_rate_quartiles_cm_per_mb"
        ],
        "physical_bp_by_distance_stratum": resource[
            "physical_bp_by_distance_stratum"
        ],
        "physical_bp_by_recombination_stratum": resource[
            "physical_bp_by_recombination_stratum"
        ],
        "physical_bp_by_joint_stratum": resource[
            "physical_bp_by_joint_stratum"
        ],
        "distance_strata": {
            key: stratum_summary(value)
            for key, value in resource["distance_strata"].items()
        },
        "recombination_strata": {
            key: stratum_summary(value)
            for key, value in resource["recombination_strata"].items()
        },
        "windows_100kb": {
            "total": len(windows),
            "nonempty": int(nonempty.sum()),
            "diversity_per_bp_quantiles": dict(
                zip(
                    ("0", "0.25", "0.5", "0.75", "1"),
                    np.quantile(diversity[nonempty], (0, 0.25, 0.5, 0.75, 1)).tolist(),
                )
            ),
            "diversity_coding_fraction_correlation": float(
                np.corrcoef(diversity[nonempty], coding[nonempty])[0, 1]
            ),
            "diversity_recombination_correlation": float(
                np.corrcoef(diversity[nonempty], recombination[nonempty])[0, 1]
            ),
            "by_contig": by_contig,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
