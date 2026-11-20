"""Prepare hap-IBD input without manufacturing haplotype phase.

Missing calls fail by default. The historical major-allele replacement remains
available explicitly for reproducing old workflows, with replacement counts.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import tempfile
from collections import Counter
from pathlib import Path


def _open_text(path: Path, mode: str):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    return opener(path, mode, encoding="utf-8")


def prepare_vcf(
    source_path: Path, output_path: Path, *, missing_policy: str = "error"
) -> dict:
    """Validate phase and publish output only after every retained record passes."""
    source_path, output_path = Path(source_path), Path(output_path)
    if missing_policy not in {"error", "drop-site", "major"}:
        raise ValueError("missing_policy must be error, drop-site or major")
    if source_path.resolve() == output_path.resolve():
        raise ValueError("Input and output must be different files")
    report = {
        "input": str(source_path),
        "output": str(output_path),
        "missing_policy": missing_policy,
        "records_read": 0,
        "records_written": 0,
        "multiallelic_records_skipped": 0,
        "nonvariant_records_skipped": 0,
        "all_missing_records_skipped": 0,
        "missing_records_dropped": 0,
        "imputed_alleles": 0,
        "unphased_homozygotes_normalized": 0,
        "scope": "Input contract only; does not infer or validate biological phase",
        "phase_set_policy": "One explicit or implicit phase set per sample/chromosome across retained heterozygous calls",
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".vcf.gz" if output_path.suffix == ".gz" else ".vcf"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=suffix, dir=output_path.parent
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    sample_names = None
    phase_sets: dict[tuple[str, str], int | None] = {}
    try:
        with (
            _open_text(source_path, "rt") as source,
            _open_text(temporary_path, "wt") as destination,
        ):
            for number, line in enumerate(source, 1):
                if line.startswith("##"):
                    if sample_names is not None:
                        raise ValueError(f"Line {number}: metadata after sample header")
                    destination.write(line)
                    continue
                if line.startswith("#CHROM\t"):
                    header = line.rstrip("\r\n").split("\t")
                    if sample_names is not None or header[:9] != [
                        "#CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                    ]:
                        raise ValueError(
                            f"Line {number}: invalid or duplicate VCF header"
                        )
                    sample_names = header[9:]
                    if (
                        not sample_names
                        or len(sample_names) != len(set(sample_names))
                        or any(not name.strip() for name in sample_names)
                    ):
                        raise ValueError("Expected nonempty, unique VCF sample IDs")
                    report["samples"] = len(sample_names)
                    destination.write(line)
                    continue
                if sample_names is None:
                    raise ValueError(f"Line {number}: missing VCF sample header")
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9 + len(sample_names):
                    raise ValueError(f"Line {number}: inconsistent VCF field count")
                report["records_read"] += 1
                if "," in fields[4]:
                    report["multiallelic_records_skipped"] += 1
                    continue
                if fields[4] == ".":
                    report["nonvariant_records_skipped"] += 1
                    continue
                formats = fields[8].split(":")
                if formats.count("GT") != 1 or len(formats) != len(set(formats)):
                    raise ValueError(
                        f"Line {number}: expected a unique GT FORMAT field"
                    )
                gt_index = formats.index("GT")
                ps_index = formats.index("PS") if "PS" in formats else None
                calls = []
                for sample_name, sample in zip(sample_names, fields[9:]):
                    parts = sample.split(":")
                    if gt_index >= len(parts) or len(parts) > len(formats):
                        raise ValueError(
                            f"Line {number}, {sample_name}: missing GT field"
                        )
                    genotype = parts[gt_index]
                    separator = "|" if "|" in genotype else "/"
                    alleles = genotype.split(separator)
                    if len(alleles) != 2 or any(
                        a not in {"0", "1", "."} for a in alleles
                    ):
                        raise ValueError(
                            f"Line {number}, {sample_name}: expected a biallelic diploid GT"
                        )
                    if separator == "/" and set(alleles) == {"0", "1"}:
                        raise ValueError(
                            f"Line {number}, {sample_name}: unphased heterozygote; "
                            "hap-IBD requires genuinely phased genotypes"
                        )
                    if "." in alleles and missing_policy == "error":
                        raise ValueError(
                            f"Line {number}, {sample_name}: missing allele; "
                            "choose an explicit --missing-policy after validating its effects"
                        )
                    calls.append((parts, separator, alleles))
                if missing_policy == "drop-site" and any(
                    "." in alleles for _, _, alleles in calls
                ):
                    report["missing_records_dropped"] += 1
                    report["all_missing_records_skipped"] += int(
                        all(a == "." for _, _, alleles in calls for a in alleles)
                    )
                    continue
                counts = Counter(
                    a for _, _, alleles in calls for a in alleles if a != "."
                )
                if not counts:
                    report["all_missing_records_skipped"] += 1
                    continue
                major = max(counts, key=lambda allele: (counts[allele], -int(allele)))
                for index, (parts, separator, alleles) in enumerate(calls, 9):
                    completed = [
                        major if allele == "." else allele for allele in alleles
                    ]
                    if separator == "/" and completed[0] != completed[1]:
                        raise ValueError(
                            f"Line {number}, {sample_names[index - 9]}: imputation cannot "
                            "establish phase for an unphased heterozygote"
                        )
                    if completed[0] != completed[1]:
                        # Homozygotes have no orientation to constrain. Sites dropped
                        # above also cannot connect phase blocks in the emitted VCF.
                        ps = (
                            parts[ps_index]
                            if ps_index is not None and ps_index < len(parts)
                            else "."
                        )
                        phase_set = None
                        if ps != ".":
                            if (
                                not ps.isascii()
                                or not ps.isdecimal()
                                or int(ps) > 2**31 - 1
                            ):
                                raise ValueError(
                                    f"Line {number}, {sample_names[index - 9]}: "
                                    "PS must be a nonnegative 32-bit integer or missing"
                                )
                            phase_set = int(ps)
                        phase_key = (fields[0], sample_names[index - 9])
                        if (
                            phase_key in phase_sets
                            and phase_sets[phase_key] != phase_set
                        ):
                            raise ValueError(
                                f"Line {number}, {sample_names[index - 9]}: unsupported "
                                f"mixed or multiple phase sets on chromosome {fields[0]}; "
                                "hap-IBD requires consistent chromosome-wide phase and "
                                "cannot connect independent PS blocks"
                            )
                        phase_sets[phase_key] = phase_set
                    report["imputed_alleles"] += alleles.count(".")
                    report["unphased_homozygotes_normalized"] += int(separator == "/")
                    parts[gt_index] = "|".join(completed)
                    fields[index] = ":".join(parts)
                destination.write("\t".join(fields) + "\n")
                report["records_written"] += 1
            if sample_names is None or report["records_written"] == 0:
                raise ValueError("No usable genotype records for hap-IBD")
        os.replace(temporary_path, output_path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--missing-policy", choices=("error", "drop-site", "major"), default="error"
    )
    args = parser.parse_args()
    print(
        json.dumps(
            prepare_vcf(args.input, args.output, missing_policy=args.missing_policy),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
