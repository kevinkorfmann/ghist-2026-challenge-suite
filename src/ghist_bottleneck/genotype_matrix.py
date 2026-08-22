from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pysam


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict_positions(positions: np.ndarray) -> tuple[np.ndarray, int, int]:
    """Give repeated integer coordinates stable sub-base ordering."""
    positions = np.asarray(positions, dtype=np.float64)
    if np.any(np.diff(positions) < 0):
        raise ValueError("positions must be nondecreasing")
    adjusted = positions.copy()
    boundaries = np.r_[0, np.flatnonzero(np.diff(positions) != 0) + 1, len(positions)]
    duplicate_records = duplicate_groups = 0
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        count = int(end - start)
        if count > 1:
            duplicate_groups += 1
            duplicate_records += count
            adjusted[start:end] += np.linspace(-0.25, 0.25, count)
    if np.any(np.diff(adjusted) <= 0):
        raise RuntimeError("sub-base ordering did not produce strictly increasing positions")
    return adjusted, duplicate_groups, duplicate_records


def export_genotype_matrix(
    vcf_path: Path,
    output_dir: Path,
    *,
    expected_records: int,
) -> dict:
    """Export every phased site as a derived-allele haplotype matrix."""
    vcf_path = vcf_path.resolve()
    output_dir = output_dir.resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    if output_dir.exists():
        raise FileExistsError(output_dir)

    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        contigs = {
            name: vcf.header.contigs[name].length for name in vcf.header.contigs
        }
    if len(contigs) != 1:
        raise ValueError(f"expected one contig, found {len(contigs)}")
    n_haplotypes = 2 * len(sample_names)
    polarization = collections.Counter()
    allele_count_histogram = collections.Counter()

    with tempfile.TemporaryDirectory(dir=output_dir.parent) as temporary_directory:
        temporary = Path(temporary_directory)
        stage = temporary / "fastcxt_input"
        stage.mkdir()
        matrix_path = stage / "genotypes.npy"
        positions_path = stage / "positions.npy"
        ordered_positions_path = stage / "strict_positions.npy"
        matrix = np.lib.format.open_memmap(
            matrix_path,
            mode="w+",
            dtype=np.uint8,
            shape=(n_haplotypes, expected_records),
        )
        positions = np.lib.format.open_memmap(
            positions_path, mode="w+", dtype=np.float64, shape=(expected_records,)
        )

        with pysam.VariantFile(str(vcf_path)) as vcf:
            for site_index, rec in enumerate(vcf):
                if site_index >= expected_records:
                    raise ValueError("VCF has more records than expected")
                if len(rec.alleles or ()) != 2 or any(len(a) != 1 for a in rec.alleles):
                    raise ValueError(f"record {site_index} is not a biallelic SNP")
                aa = rec.info.get("AA")
                if isinstance(aa, tuple):
                    aa = aa[0]
                if aa == rec.ref:
                    flip = False
                    polarization["ancestral_is_reference"] += 1
                elif aa == rec.alts[0]:
                    flip = True
                    polarization["ancestral_is_alternate"] += 1
                else:
                    raise ValueError(f"record {site_index} has invalid ancestral allele {aa!r}")

                column = np.empty(n_haplotypes, dtype=np.uint8)
                for individual_index, call in enumerate(rec.samples.values()):
                    gt = call.get("GT")
                    if gt is None or len(gt) != 2 or any(a is None or a < 0 for a in gt):
                        raise ValueError(f"record {site_index} has missing or non-diploid GT")
                    if not call.phased:
                        raise ValueError(f"record {site_index} has an unphased GT")
                    column[2 * individual_index : 2 * individual_index + 2] = gt
                if flip:
                    column = 1 - column
                derived_count = int(column.sum())
                if not 0 < derived_count < n_haplotypes:
                    raise ValueError(f"record {site_index} is fixed after polarization")
                matrix[:, site_index] = column
                positions[site_index] = rec.pos
                allele_count_histogram[derived_count] += 1
            records = site_index + 1
        if records != expected_records:
            raise ValueError(f"VCF has {records} records, expected {expected_records}")
        matrix.flush()
        positions.flush()

        adjusted, duplicate_groups, duplicate_records = strict_positions(positions)
        np.save(ordered_positions_path, adjusted)
        metadata = {
            "source": str(vcf_path),
            "source_sha256": sha256(vcf_path),
            "records": records,
            "individuals": len(sample_names),
            "haplotypes": n_haplotypes,
            "sample_names": sample_names,
            "contigs": contigs,
            "genotypes_shape": [n_haplotypes, records],
            "genotypes_dtype": "uint8",
            "polarization": dict(polarization),
            "derived_allele_count_histogram": {
                str(k): v for k, v in sorted(allele_count_histogram.items())
            },
            "duplicate_position_groups": duplicate_groups,
            "records_at_duplicate_positions": duplicate_records,
            "strict_position_rule": "stable offsets from -0.25 to 0.25 bp within each repeated coordinate",
        }
        (stage / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        shutil.move(stage, output_dir)

    metadata["files"] = {
        name: {
            "bytes": (output_dir / name).stat().st_size,
            "sha256": sha256(output_dir / name),
        }
        for name in ("genotypes.npy", "positions.npy", "strict_positions.npy")
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-records", type=int, required=True)
    args = parser.parse_args()
    result = export_genotype_matrix(
        args.vcf, args.output_dir, expected_records=args.expected_records
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
