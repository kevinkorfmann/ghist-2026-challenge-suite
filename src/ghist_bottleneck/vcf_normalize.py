from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import pysam


AA_HEADER = (
    '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. '
    'Format: AA. AA: Ancestral allele">\n'
)


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def body_digest(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    records = 0
    with gzip.open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b"#"):
                continue
            digest.update(line)
            records += 1
    return digest.hexdigest(), records


def normalize_vcf(source: Path, destination: Path, report_path: Path) -> dict:
    """Repair the known malformed metadata line without changing any VCF record."""
    source = source.resolve()
    destination = destination.resolve()
    report_path = report_path.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffixes[-2:] != [".vcf", ".gz"]:
        raise ValueError("destination must end in .vcf.gz")

    source_body_sha256, source_records = body_digest(source)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary_directory:
        temporary = Path(temporary_directory)
        plain = temporary / destination.with_suffix("").name
        repaired = 0
        column_header_seen = False
        with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle:
            with plain.open("wt", encoding="utf-8", newline="") as output_handle:
                for line in input_handle:
                    if line.startswith("##INFO=<ID=AA,"):
                        output_handle.write(AA_HEADER)
                        repaired += 1
                    else:
                        output_handle.write(line)
                    column_header_seen |= line.startswith("#CHROM\t")
        if repaired != 1:
            raise ValueError(f"expected exactly one AA metadata line, found {repaired}")
        if not column_header_seen:
            raise ValueError("VCF column header is absent")

        compressed_path = Path(
            pysam.tabix_index(str(plain), preset="vcf", force=True, keep_original=False)
        )
        index_path = Path(str(compressed_path) + ".tbi")
        normalized_body_sha256, normalized_records = body_digest(compressed_path)
        if source_records != normalized_records or source_body_sha256 != normalized_body_sha256:
            raise RuntimeError("normalization changed one or more VCF records")

        with pysam.VariantFile(str(compressed_path)) as normalized:
            samples = list(normalized.header.samples)
            contigs = {
                name: normalized.header.contigs[name].length
                for name in normalized.header.contigs
            }

        shutil.move(compressed_path, destination)
        shutil.move(index_path, Path(str(destination) + ".tbi"))

    report = {
        "source": str(source),
        "destination": str(destination),
        "repair": "canonicalized malformed AA metadata line",
        "repaired_header_lines": repaired,
        "records": source_records,
        "record_body_sha256": source_body_sha256,
        "samples": len(samples),
        "contigs": contigs,
        "source_bytes": source.stat().st_size,
        "source_md5": file_digest(source, "md5"),
        "source_sha256": file_digest(source, "sha256"),
        "destination_bytes": destination.stat().st_size,
        "destination_md5": file_digest(destination, "md5"),
        "destination_sha256": file_digest(destination, "sha256"),
        "index_bytes": Path(str(destination) + ".tbi").stat().st_size,
        "index_sha256": file_digest(Path(str(destination) + ".tbi"), "sha256"),
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(normalize_vcf(args.source, args.destination, args.report), indent=2))


if __name__ == "__main__":
    main()
