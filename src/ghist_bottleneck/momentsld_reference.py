from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile

import joblib
import numpy as np

from .features import LD_EDGES
from .sfs_baseline import projected_sfs_from_vcf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_hdf5(vcf: str | Path) -> Path:
    import allel

    vcf = Path(vcf)
    hdf5 = Path(str(vcf).split(".vcf")[0] + ".h5")
    allel.vcf_to_hdf5(
        str(vcf),
        str(hdf5),
        fields="*",
        exclude_fields=["calldata/GQ"],
        overwrite=True,
    )
    return hdf5


def compute_region(
    vcf: str | Path,
    out: str | Path,
    *,
    block_index: int,
    block_length: int = 2_000_000,
) -> dict:
    if not 0 <= block_index < 50:
        raise ValueError("block_index must be in [0, 50)")
    import moments.LD as LD

    left = block_index * block_length
    right = left + block_length
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed") as bed:
        bed.write(f"1\t{left}\t{right}\n")
        bed.flush()
        statistics = LD.Parsing.compute_ld_statistics(
            str(vcf),
            bed_file=bed.name,
            chromosome="1",
            bp_bins=LD_EDGES.tolist(),
            use_genotypes=False,
            use_h5=True,
            report=True,
            use_cache=True,
        )
    record = {
        "block_index": block_index,
        "block_start": left,
        "block_end": right,
        "statistics": statistics,
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(record, out)
    return record


def fit_regions(
    region_dir: str | Path,
    vcf: str | Path,
    *,
    restarts: int = 32,
    seed: int = 7411,
    mutation_rate: float = 3.42e-8,
    recombination_rate: float = 1.87e-8,
    sequence_length: float = 100_000_000,
    joint_sfs: bool = True,
) -> dict:
    import moments
    import moments.LD as LD

    region_dir = Path(region_dir)
    paths = sorted(region_dir.glob("region_*.joblib"))
    if len(paths) != 50:
        raise ValueError(f"expected 50 LD regions, found {len(paths)}")
    records = [joblib.load(path) for path in paths]
    indices = [int(record["block_index"]) for record in records]
    if sorted(indices) != list(range(50)) or len(set(indices)) != 50:
        raise ValueError("LD region indices must be exactly 0 through 49")
    all_data = {index: record["statistics"] for index, record in zip(indices, records)}
    np.random.seed(seed)
    bootstrapped = LD.Parsing.bootstrap_data(all_data)
    raw_edges, means, varcovs, statistics = LD.Parsing.subset_data(
        bootstrapped, ["ALL"]
    )
    recombination_edges = raw_edges * recombination_rate
    model_functions = [LD.Demographics1D.two_epoch]
    data = [means, varcovs]
    spectrum_metadata = None
    if joint_sfs:
        spectrum, spectrum_metadata = projected_sfs_from_vcf(vcf, polarized=True)
        data.append(moments.Spectrum(spectrum))
        model_functions.append(moments.Demographics1D.two_epoch)

    lower = np.array([0.02, 1e-4, 5_000.0])
    upper = np.array([2.0, 5.0, 50_000.0])
    rng = np.random.default_rng(seed)
    fits = []
    failures = []
    for restart in range(restarts):
        start = np.exp(rng.uniform(np.log(lower), np.log(upper)))
        try:
            estimate, objective = LD.Inference.optimize_log_lbfgsb(
                start,
                data,
                model_functions,
                rs=recombination_edges,
                u=mutation_rate,
                lower_bound=lower,
                upper_bound=upper,
                use_afs=joint_sfs,
                Leff=sequence_length if joint_sfs else None,
                multinom=False,
                ns=[114] if joint_sfs else None,
                statistics=statistics,
                maxiter=1000,
            )
            estimate = np.asarray(estimate, dtype=float)
            if np.isfinite([*estimate, objective]).all():
                fits.append((float(objective), estimate, restart))
            else:
                failures.append({"restart": restart, "reason": "nonfinite result"})
        except Exception as error:
            failures.append(
                {"restart": restart, "reason": f"{type(error).__name__}: {error}"}
            )
    if not fits:
        raise RuntimeError(f"all momentsLD restarts failed: {failures[:3]}")
    objective, estimate, best_restart = min(fits, key=lambda item: item[0])
    fraction, scaled_time, ancestral_ne = estimate
    generations = scaled_time * 2 * ancestral_ne
    return {
        "generations": float(generations),
        "post_decline_fraction": float(fraction),
        "ancestral_effective_size": float(ancestral_ne),
        "scaled_time_2N": float(scaled_time),
        "negative_composite_log_likelihood": objective,
        "joint_sfs": joint_sfs,
        "best_restart": int(best_restart),
        "successful_restarts": len(fits),
        "failed_restarts": failures,
        "restarts": restarts,
        "bootstrap_regions": 50,
        "bootstrap_seed": seed,
        "bp_bin_edges": raw_edges.tolist(),
        "recombination_rate": recombination_rate,
        "mutation_rate": mutation_rate,
        "sequence_length": sequence_length,
        "sfs": spectrum_metadata,
        "region_sha256": {path.name: sha256(path) for path in paths},
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--vcf", type=Path, required=True)
    region = subparsers.add_parser("region")
    region.add_argument("--vcf", type=Path, required=True)
    region.add_argument("--out", type=Path, required=True)
    region.add_argument("--block-index", type=int, required=True)
    fit = subparsers.add_parser("fit")
    fit.add_argument("--vcf", type=Path, required=True)
    fit.add_argument("--region-dir", type=Path, required=True)
    fit.add_argument("--out", type=Path, required=True)
    fit.add_argument("--restarts", type=int, default=32)
    fit.add_argument("--seed", type=int, default=7411)
    fit.add_argument("--ld-only", action="store_true")
    args = parser.parse_args()
    if args.command == "prepare":
        path = prepare_hdf5(args.vcf)
        print(json.dumps({"hdf5": str(path), "sha256": sha256(path)}, indent=2))
    elif args.command == "region":
        result = compute_region(args.vcf, args.out, block_index=args.block_index)
        print(json.dumps({key: result[key] for key in result if key != "statistics"}, indent=2))
    else:
        result = fit_regions(
            args.region_dir,
            args.vcf,
            restarts=args.restarts,
            seed=args.seed,
            joint_sfs=not args.ld_only,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
