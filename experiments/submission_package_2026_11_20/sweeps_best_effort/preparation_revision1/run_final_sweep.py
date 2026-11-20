"""Frozen best-effort deadline inference; no calibration, truth access, or submission."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pysam

from ghist_single_sweep.background import callable_raisd_windows, masked_variants
from ghist_single_sweep.peaks import read_raisd
from ghist_single_sweep.predict import predict_single_sweep
from ghist_single_sweep.scan import (
    merge_intervals,
    read_bed_intervals,
    read_phased_biallelic_vcf,
)
from ghist_single_sweep.soft import snp_h12_peak
from ghist_submission.validate import validate_submission
from predict_task7_sweep_bgs import run as run_task7
from predict_task7_sweep_bgs import validate_raisd_site_count, write_masked_vcf

PROTOCOL_SHA256 = "acc0ce4de87e4e99c282b97947dbe87b5748a2d8034dd519eb33410237a4fd17"
SOURCE_PATHS = [
    "scripts/predict_task7_sweep_bgs.py",
    "src/ghist_single_sweep/background.py",
    "src/ghist_single_sweep/scan.py",
    "src/ghist_single_sweep/peaks.py",
    "src/ghist_single_sweep/predict.py",
    "src/ghist_single_sweep/recovery.py",
    "src/ghist_single_sweep/soft.py",
    "src/ghist_submission/validate.py",
    "src/ghist_single_sweep/__init__.py",
    "src/ghist_submission/__init__.py",
    "src/ghist_relatedness/__init__.py",
    "src/ghist_relatedness/labels.py",
]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_manifest(root):
    return {p: sha256(root / p) for p in SOURCE_PATHS}


def intersection(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = (
            max(int(left[i][0]), int(right[j][0])),
            min(int(left[i][1]), int(right[j][1])),
        )
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return np.asarray(out, dtype=np.int64).reshape(-1, 2)


def nearest(position, intervals):
    if not len(intervals):
        raise ValueError("No eligible intervals for nearest callable base")
    candidates = [min(max(int(position), int(a)), int(b) - 1) for a, b in intervals]
    return min(candidates, key=lambda x: (abs(x - position), x))


def final_geometry(centers, half_width, mask, length, multiple=False):
    records = []
    intervals = []
    for center in centers:
        center = int(center)
        adjusted = nearest(center, mask)
        component = next((int(a), int(b)) for a, b in mask if a <= adjusted < b)
        interval = [
            max(0, component[0], adjusted - half_width),
            min(length, component[1], adjusted + half_width),
        ]
        if not 0 <= interval[0] <= adjusted < interval[1] <= length:
            raise ValueError("Invalid final callable containment")
        intervals.append(interval)
        records.append(
            {
                "raw_center": center,
                "final_center": adjusted,
                "callable_component": component,
                "interval": interval,
                "center_adjusted": adjusted != center,
                "boundary_clipped": interval
                != [adjusted - half_width, adjusted + half_width],
            }
        )
    if not intervals:
        raise ValueError("No final intervals")
    if multiple:
        intervals = merge_intervals(np.asarray(intervals)).tolist()
    elif len(intervals) != 1:
        raise ValueError("Single-sweep task requires one interval")
    return intervals, records


def task8_candidates(report, reference, exons):
    """A declared transferred-score heuristic, never a calibrated null."""
    order = np.lexsort((report[:, 0], -report[:, 6]))
    reference_order = (
        np.lexsort((reference[:, 0], -reference[:, 6])) if len(reference) else []
    )
    background = None
    if len(reference_order):
        anchor = reference[reference_order[0], 0] - 1
        distant = reference[np.abs(reference[:, 0] - 1 - anchor) > 2_000_000]
        if len(distant):
            background = float(np.max(distant[:, 6]))
    threshold = None if background is None else 2 * background
    chosen = []
    selected_positions = []
    if threshold is not None:
        for index in order:
            position, score = int(report[index, 0]) - 1, float(report[index, 6])
            if score <= threshold:
                continue
            if any(
                abs(position - item["raisd_position"]) < 1_000_000 for item in chosen
            ):
                continue
            selected_positions.append(position)
            snapped = nearest(position, exons) if len(exons) else position
            if not len(exons) or abs(snapped - position) > 60_000:
                continue
            chosen.append(
                {
                    "raisd_position": position,
                    "center": snapped,
                    "score": score,
                    "exon_snap": snapped - position,
                }
            )
    reason = None
    if not chosen:
        reason = (
            "missing_reference_background"
            if background is None
            else "no_above_threshold_callable_exonic_candidates"
        )
        index = order[0]
        position = int(report[index, 0]) - 1
        snapped = nearest(position, exons) if len(exons) else position
        if abs(snapped - position) > 60_000:
            snapped = position
        chosen = [
            {
                "raisd_position": position,
                "center": snapped,
                "score": float(report[index, 6]),
                "exon_snap": snapped - position,
            }
        ]
    return {
        "method": "provisional_legacy_raisd_transfer_v1",
        "candidates": chosen,
        "reference_background_score": background,
        "threshold": threshold,
        "fallback_reason": reason,
        "calibrated_null": False,
        "scientifically_validated": False,
        "truth_used_for_inference": False,
    }


def input_contract(vcf, mask_path):
    with pysam.VariantFile(str(vcf)) as handle:
        contigs = list(handle.header.contigs)
        if len(contigs) != 1 or contigs[0].removeprefix("chr") != "9":
            raise ValueError("Final sweep input must contain exactly chromosome9")
        contig = contigs[0]
        length = int(handle.header.contigs[contig].length)
        samples = len(handle.header.samples)
    mask = read_bed_intervals(mask_path, contig)
    if length <= 0 or mask[-1, 1] > length:
        raise ValueError("Final accessibility outside exact declared contig")
    return contig, length, samples, mask


def checked_reference(path, source, mask_hash, length):
    completion_path = path / "completion.json"
    completion = json.loads(completion_path.read_text())
    if completion.get("state") != "COMPLETED" or completion.get("task") != 7:
        raise ValueError("Task8 requires completed Task7 inference reference")
    if (
        completion.get("protocol_sha256") != PROTOCOL_SHA256
        or completion.get("source_sha256") != source
    ):
        raise ValueError("Task7 reference method/source mismatch")
    if (
        completion["inputs"]["mask"]["sha256"] != mask_hash
        or completion["contig_length"] != length
    ):
        raise ValueError("Task7 reference callable exposure differs")
    for relative, digest in completion["artifact_sha256"].items():
        if sha256(path / relative) != digest:
            raise ValueError("Task7 reference artifact changed")
    return path / "raw" / "RAiSD_Report.task7", sha256(completion_path)


def run(args):
    protocol = Path(__file__).with_name("method_protocol.json")
    freeze = json.loads(args.source_freeze.read_text())
    if (
        sha256(protocol) != PROTOCOL_SHA256
        or freeze["runner_sha256"] != sha256(__file__)
        or freeze["protocol_sha256"] != PROTOCOL_SHA256
    ):
        raise ValueError("Runner/protocol freeze mismatch")
    for module_name, relative in {
        "ghist_single_sweep.scan": "src/ghist_single_sweep/scan.py",
        "ghist_single_sweep.predict": "src/ghist_single_sweep/predict.py",
        "ghist_single_sweep.soft": "src/ghist_single_sweep/soft.py",
        "ghist_single_sweep.background": "src/ghist_single_sweep/background.py",
        "ghist_submission.validate": "src/ghist_submission/validate.py",
        "predict_task7_sweep_bgs": "scripts/predict_task7_sweep_bgs.py",
    }.items():
        if (
            Path(importlib.import_module(module_name).__file__).resolve()
            != (args.source_root / relative).resolve()
        ):
            raise ValueError("Imported source is outside the frozen source root")
    source = source_manifest(args.source_root)
    if source != freeze["source_sha256"]:
        raise ValueError("Source differs from pre-input frozen methods")
    input_paths = {"vcf": args.vcf, "mask": args.accessible}
    expected = {"vcf": args.vcf_sha256, "mask": args.mask_sha256}
    if args.task == 8:
        if (
            args.exons is None
            or args.exons_sha256 is None
            or args.task7_reference is None
        ):
            raise ValueError("Task8 needs final exons and completed Task7 reference")
        input_paths["exons"] = args.exons
        expected["exons"] = args.exons_sha256
    inputs = {
        name: {"path": str(path.resolve()), "sha256": sha256(path)}
        for name, path in input_paths.items()
    }
    if any(inputs[name]["sha256"] != expected[name] for name in inputs):
        raise ValueError("Final input differs from authenticated download hash")
    if args.out_dir.exists():
        raise ValueError("Refusing to overwrite deadline inference artifacts")
    args.out_dir.mkdir(parents=True)
    started = time.time()
    runtime = {
        "python": platform.python_version(),
        "executable_sha256": sha256(sys.executable),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "pysam", "PyYAML")
        },
        "raisd_binary_sha256": sha256(args.raisd),
    }
    bound = {
        "task": args.task,
        "protocol_sha256": PROTOCOL_SHA256,
        "source_freeze_sha256": sha256(args.source_freeze),
        "runner_sha256": sha256(__file__),
        "source_sha256": source,
        "inputs": inputs,
        "runtime": runtime,
        "started_unix": started,
        "scientifically_validated": False,
        "release_ready": False,
        "external_submissions": 0,
    }
    dump(args.out_dir / "started.json", bound)
    try:
        contig, length, samples, mask = input_contract(args.vcf, args.accessible)
        raw = args.out_dir / "raw"
        reference_hash = None
        if args.task == 7:
            result = run_task7(args.vcf, args.accessible, "final", args.raisd, raw)
            centers = [result["center"]]
        else:
            raw.mkdir()
            variants = read_phased_biallelic_vcf(args.vcf)
            if args.task == 5:
                result = predict_single_sweep(variants, mask)
                centers = [result["center"]]
            elif args.task == 6:
                result = snp_h12_peak(variants, mask, 50, 0)
                result.update(method="h12_snps50_maf0", scientifically_validated=False)
                centers = [result["position"]]
            else:
                variants, mask = masked_variants(variants, mask)
                if len(variants.positions) < 50:
                    raise ValueError("Fewer than50callable SNPs")
                canonical = raw / "masked.vcf"
                write_masked_vcf(variants, canonical)
                with (raw / "raisd.log").open("w") as log:
                    subprocess.run(
                        [
                            str(args.raisd.resolve()),
                            "-n",
                            "task8",
                            "-I",
                            str(canonical.resolve()),
                            "-w",
                            "50",
                            "-R",
                        ],
                        cwd=raw,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                validate_raisd_site_count(
                    raw / "RAiSD_Info.task8", len(variants.positions)
                )
                report = callable_raisd_windows(
                    read_raisd(raw / "RAiSD_Report.task8"), mask, length
                )
                reference_path, reference_hash = checked_reference(
                    args.task7_reference, source, inputs["mask"]["sha256"], length
                )
                reference = callable_raisd_windows(
                    read_raisd(reference_path), mask, length
                )
                exons = read_bed_intervals(args.exons, contig)
                if exons[-1, 1] > length:
                    raise ValueError("Final exon annotation outside contig")
                exons = intersection(exons, mask)
                result = task8_candidates(report, reference, exons)
                centers = [row["center"] for row in result["candidates"]]
                result["task7_completion_sha256"] = reference_hash
                result["callable_exon_intervals"] = len(exons)
            dump(raw / "prediction.json", result)
            del variants
        half_width = {5: 60000, 6: 75000, 7: 80000, 8: 60000}[args.task]
        intervals, adjustments = final_geometry(
            centers, half_width, mask, length, args.task == 8
        )
        bed = args.out_dir / "deadline_submission.bed"
        bed.write_text("".join(f"9\t{a}\t{b}\n" for a, b in intervals))
        validation = validate_submission(bed, args.task, "final")
        if (
            any(
                sha256(path) != inputs[name]["sha256"]
                for name, path in input_paths.items()
            )
            or source_manifest(args.source_root) != source
            or sha256(__file__) != bound["runner_sha256"]
            or sha256(protocol) != PROTOCOL_SHA256
            or sha256(args.source_freeze) != bound["source_freeze_sha256"]
            or sha256(args.raisd) != runtime["raisd_binary_sha256"]
        ):
            raise ValueError("Input/source/runtime file changed during inference")
        if (
            reference_hash is not None
            and checked_reference(
                args.task7_reference, source, inputs["mask"]["sha256"], length
            )[1]
            != reference_hash
        ):
            raise ValueError("Task7 reference changed during Task8 inference")
        artifacts = {
            str(path.relative_to(args.out_dir)): sha256(path)
            for path in sorted(args.out_dir.rglob("*"))
            if path.is_file()
        }
        completion = {
            **bound,
            "state": "COMPLETED",
            "completed_unix": time.time(),
            "duration_seconds": time.time() - started,
            "contig": contig,
            "contig_length": length,
            "samples_diploid": samples,
            "final_intervals": intervals,
            "geometry_adjustments": adjustments,
            "format_validation": validation,
            "artifact_sha256": artifacts,
            "limitation": "Best-effort deadline candidate only. Scientific qualification and matched calibration absent; original validation failures remain.",
        }
        dump(args.out_dir / "completion.json", completion)
        return completion
    except BaseException as error:
        dump(
            args.out_dir / "failure.json",
            {
                **bound,
                "state": "FAILED",
                "error_type": type(error).__name__,
                "error": str(error),
                "ended_unix": time.time(),
            },
        )
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", type=int, choices=(5, 6, 7, 8), required=True)
    for name in (
        "vcf",
        "accessible",
        "raisd",
        "source-root",
        "source-freeze",
        "out-dir",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("vcf-sha256", "mask-sha256"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--exons", type=Path)
    parser.add_argument("--exons-sha256")
    parser.add_argument("--task7-reference", type=Path)
    result = run(parser.parse_args())
    print(
        json.dumps(
            {
                "task": result["task"],
                "state": result["state"],
                "duration_seconds": result["duration_seconds"],
                "scientifically_validated": False,
            }
        )
    )


if __name__ == "__main__":
    main()
