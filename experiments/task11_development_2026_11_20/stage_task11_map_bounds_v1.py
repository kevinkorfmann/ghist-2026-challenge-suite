"""Deterministic VCF staging to exact map bounds; retained bytes are unchanged."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bounds_from_map(path):
    rows = {}
    with Path(path).open() as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 4:
                raise ValueError("Expected four PLINK map columns")
            chromosome, _, cm_text, bp_text = fields
            bp, cm = float(bp_text), float(cm_text)
            if not math.isfinite(bp) or bp <= 0 or not bp.is_integer() or not math.isfinite(cm) or cm < 0:
                raise ValueError("Invalid map position")
            bp = int(bp)
            if chromosome in rows and (bp <= rows[chromosome][-1][0] or cm < rows[chromosome][-1][1]):
                raise ValueError("Map positions must increase and cM must not decrease")
            rows.setdefault(chromosome, []).append((bp, cm))
    if not rows or any(len(values) < 2 for values in rows.values()):
        raise ValueError("Each map chromosome needs at least two positions")
    return {chromosome: (values[0][0], values[-1][0]) for chromosome, values in rows.items()}


def binary_reader(path):
    with Path(path).open("rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if compressed else Path(path).open("rb")


def stage(source, map_path, target, expected_map_sha, min_markers=50):
    source, map_path, target = Path(source), Path(map_path), Path(target)
    if not isinstance(min_markers, int) or min_markers < 1:
        raise ValueError("Minimum raw retained marker count must be positive")
    if sha(map_path) != expected_map_sha:
        raise ValueError("Frozen map hash mismatch")
    input_sha = sha(source)
    bounds = bounds_from_map(map_path)
    counts = {chromosome: {"map_min_bp": low, "map_max_bp": high, "input_records": 0,
        "retained_raw_markers": 0, "excluded_below": 0, "excluded_above": 0,
        "input_min_bp": None, "input_max_bp": None, "retained_min_bp": None, "retained_max_bp": None}
        for chromosome, (low, high) in bounds.items()}
    input_lines, retained_lines, all_output = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
    samples, previous, seen = None, {}, set()
    # Exclusive creation prevents overwriting an original input or prior stage.
    with (binary_reader(source) as reader, target.open("xb") as raw,
          gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as writer):
        for line in reader:
            if line.startswith(b"##") and samples is None:
                writer.write(line)
                all_output.update(line)
                continue
            if line.startswith(b"#CHROM\t") and samples is None:
                fields = line.rstrip(b"\r\n").split(b"\t")
                if fields[:9] != [b"#CHROM", b"POS", b"ID", b"REF", b"ALT", b"QUAL", b"FILTER", b"INFO", b"FORMAT"]:
                    raise ValueError("Invalid VCF columns")
                samples = fields[9:]
                if len(samples) < 8 or len(set(samples)) != len(samples):
                    raise ValueError("Expected at least eight unique sample IDs")
                writer.write(line)
                all_output.update(line)
                continue
            if samples is None or line.startswith(b"#"):
                raise ValueError("Missing or misplaced VCF header")
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 9 + len(samples):
                raise ValueError("Genotype column count differs from roster")
            chromosome, position = fields[0].decode("utf-8"), int(fields[1])
            if chromosome not in bounds:
                raise ValueError(f"VCF chromosome absent from frozen map: {chromosome}")
            if position <= 0 or position < previous.get(chromosome, 0):
                raise ValueError("VCF positions must be positive and nondecreasing within chromosome")
            previous[chromosome] = position
            seen.add(chromosome)
            entry = counts[chromosome]
            entry["input_records"] += 1
            entry["input_min_bp"] = position if entry["input_min_bp"] is None else min(entry["input_min_bp"], position)
            entry["input_max_bp"] = position if entry["input_max_bp"] is None else max(entry["input_max_bp"], position)
            input_lines.update(line)
            low, high = bounds[chromosome]
            if position < low:
                entry["excluded_below"] += 1
            elif position > high:
                entry["excluded_above"] += 1
            else:
                entry["retained_raw_markers"] += 1
                entry["retained_min_bp"] = position if entry["retained_min_bp"] is None else min(entry["retained_min_bp"], position)
                entry["retained_max_bp"] = position if entry["retained_max_bp"] is None else max(entry["retained_max_bp"], position)
                retained_lines.update(line)
                all_output.update(line)
                writer.write(line)
    if samples is None or seen != set(bounds):
        raise ValueError("VCF must cover every frozen map chromosome")
    if any(entry["retained_raw_markers"] < min_markers for entry in counts.values()):
        raise ValueError("Too few raw retained markers on a map chromosome")
    if sha(source) != input_sha or sha(map_path) != expected_map_sha:
        raise ValueError("Input/map changed during staging")
    return {"schema": "task11-map-bound-staging-v1", "source_vcf": str(source.resolve()),
        "source_vcf_sha256": input_sha, "map": str(map_path.resolve()), "map_sha256": expected_map_sha,
        "output_vcf": str(target.resolve()), "output_vcf_sha256": sha(target),
        "inclusive_map_bounds": True, "samples": len(samples), "chromosomes": counts,
        "input_records": sum(v["input_records"] for v in counts.values()),
        "retained_records": sum(v["retained_raw_markers"] for v in counts.values()),
        "excluded_records": sum(v["excluded_below"] + v["excluded_above"] for v in counts.values()),
        "minimum_retained_raw_markers": min_markers,
        "input_variant_lines_sha256": input_lines.hexdigest(), "retained_variant_lines_sha256": retained_lines.hexdigest(),
        "output_uncompressed_sha256": all_output.hexdigest(), "staging_script_sha256": sha(__file__),
        "scope": "Only out-of-map variant records removed. Every retained line/header/genotype/phase/FORMAT byte preserved. No cM extrapolation or imputation. Frozen adapter applies its existing preparation afterward."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("vcf", "map", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--map-sha256", required=True)
    parser.add_argument("--min-markers", type=int, default=50)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    audit = stage(args.vcf, args.map, args.out / "map_bounded.vcf.gz", args.map_sha256, args.min_markers)
    (args.out / "map_crop_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"retained": audit["retained_records"], "excluded": audit["excluded_records"], "output_sha256": audit["output_vcf_sha256"]}))


if __name__ == "__main__":
    main()
