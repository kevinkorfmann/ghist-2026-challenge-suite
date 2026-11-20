"""BAM -> ANGSD SAF -> globally folded realSFS -> demographic fit.

No uncertain hard genotype calls, SNP ascertainment, or assumed ancestral reference.
The pipeline deliberately has no submission writer or final-data downloader.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
import subprocess

import numpy as np


TASKS = {
    "bottleneck": {"samples": {"pop1": 20}, "sample_sizes": [40], "depth": 2.0},
    "split_isolation": {"samples": {"pop1": 10, "pop2": 10}, "sample_sizes": [20, 20], "depth": 4.0},
}
MUTATION_RATE = 1e-7
SEQUENCE_LENGTH = 10_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")


def audit_inputs(directory: Path, reference: Path, task: str, expected_length: int = SEQUENCE_LENGTH) -> dict:
    """Require the full roster and a sorted, indexed, one-contig BAM for each SM."""
    import pysam

    if task not in TASKS or expected_length <= 0:
        raise ValueError("unknown task or invalid expected length")
    directory, reference = Path(directory).resolve(), Path(reference).resolve()
    samples = {}
    for line in (directory / "popfile.txt").read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 2 or fields[0] in samples:
            raise ValueError("popfile must have exactly two columns and unique samples")
        samples[fields[0]] = fields[1]
    if task == "bottleneck":
        if len(samples) != 20 or len(set(samples.values())) != 1:
            raise ValueError("bottleneck requires 20 samples from one population")
        populations = sorted(set(samples.values()))
    else:
        observed = {p: sum(v == p for v in samples.values()) for p in set(samples.values())}
        if observed != TASKS[task]["samples"]:
            raise ValueError("split requires exactly ten pop1 and ten pop2 samples")
        populations = ["pop1", "pop2"]
    with pysam.FastaFile(str(reference)) as fasta:
        if len(fasta.references) != 1 or tuple(fasta.lengths) != (expected_length,):
            raise ValueError("reference must contain one contig at the expected length")
        contig = fasta.references[0]
        # Do not treat ambiguous reference positions as callable scale.
        non_acgt = sum(sum(base not in "ACGT" for base in fasta.fetch(contig, start, min(start + 1_000_000, expected_length)).upper())
                       for start in range(0, expected_length, 1_000_000))
        if non_acgt:
            raise ValueError("reference contains ambiguous bases; explicit accessibility handling is required")
    records = {}
    for bam in sorted(directory.glob("*.bam")):
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            header = handle.header.to_dict()
            names = {rg.get("SM") for rg in header.get("RG", [])}
            if len(names) != 1 or None in names:
                raise ValueError(f"BAM requires one unambiguous SM: {bam.name}")
            sample = next(iter(names))
            if sample not in samples or sample in records:
                raise ValueError("BAM roster differs from popfile or has duplicate SM")
            if header.get("HD", {}).get("SO") != "coordinate" or not handle.has_index():
                raise ValueError("BAM must be coordinate sorted and indexed")
            if handle.references != (contig,) or handle.lengths != (expected_length,):
                raise ValueError("BAM contig dictionary differs from reference")
            records[sample] = {"path": str(bam), "sha256": sha256(bam), "population": samples[sample]}
    if set(records) != set(samples):
        raise ValueError("BAM roster is incomplete")
    return {"task": task, "reference": str(reference), "reference_sha256": sha256(reference),
            "popfile_sha256": sha256(directory / "popfile.txt"), "length": expected_length, "contig": contig,
            "population_order": populations, "sample_sizes": TASKS[task]["sample_sizes"], "samples": records}


def folded_spectrum(values: np.ndarray, sample_sizes: list[int]):
    """Import realSFS global folding, symmetrizing the half-frequency boundary.

    realSFS keeps every bin with total count <= half the pooled chromosomes.
    Its distinct equal-frequency complementary bins are unidentifiable and can
    differ numerically. moments divides their sum equally, so do that once.
    Do not call Spectrum.fold() on this already-folded output.
    """
    import moments

    shape = tuple(n + 1 for n in sample_sizes)
    data = np.asarray(values, dtype=float)
    if data.size != int(np.prod(shape)) or np.any(~np.isfinite(data)) or np.any(data < 0):
        raise ValueError("invalid realSFS count array")
    data = data.reshape(shape).copy()
    counts = np.indices(shape).sum(axis=0)
    outside = counts * 2 > sum(sample_sizes)
    if np.any(data[outside] > 1e-5):
        raise ValueError("SFS is not globally folded")
    boundary = counts * 2 == sum(sample_sizes)
    reverse = data[tuple(slice(None, None, -1) for _ in shape)]
    data[boundary] = (data[boundary] + reverse[boundary]) / 2
    mask = outside.copy()
    mask[(0,) * len(shape)] = True
    return moments.Spectrum(data, mask=mask, data_folded=True, mask_corners=True)


def expected_spectrum(parameters, sample_sizes, task):
    import moments

    if task == "bottleneck":
        return moments.Demographics1D.two_epoch(parameters, sample_sizes).fold()
    if task == "split_isolation":
        return moments.Demographics2D.split_mig([*parameters, 0], sample_sizes).fold()
    raise ValueError("unknown task")


def physical_parameters(task: str, parameters, theta: float, callable_length: float) -> dict:
    if not np.isfinite(theta) or theta <= 0 or not np.isfinite(callable_length) or callable_length <= 0:
        raise ValueError("positive finite theta and callable length required")
    ancestral = theta / (4 * MUTATION_RATE * callable_length)
    if task == "bottleneck":
        fraction, time = parameters
        target = {"generations": 2 * ancestral * time, "post_decline_fraction": fraction}
    elif task == "split_isolation":
        size1, size2, time = parameters
        target = {"pop1_population_size": ancestral * size1, "pop2_population_size": ancestral * size2,
                  "generations_since_split": 2 * ancestral * time}
    else:
        raise ValueError("unknown task")
    return {"ancestral_population_size": ancestral, "targets": target}


def fit_spectrum(values, sample_sizes, task, *, restarts: int = 8, seed: int = 2026112013) -> dict:
    """Independent fixed optimizer starts; no adjustment using challenge answers."""
    import moments
    from scipy.optimize import minimize

    if restarts < 2:
        raise ValueError("at least two independent optimization starts required")
    if list(sample_sizes) != TASKS[task]["sample_sizes"]:
        raise ValueError("sample dimensions differ from the challenge contract")
    data = folded_spectrum(values, sample_sizes)
    callable_length = float(np.sum(values))
    if float(data.sum()) <= 0 or callable_length <= 0:
        raise ValueError("no segregating information")
    bounds = [(.001, .9999), (.0001, 10.0)] if task == "bottleneck" else [(.05, 20.0), (.05, 20.0), (.0001, 10.0)]
    rng = np.random.default_rng(seed)
    starts = [np.array([.25, .1])] if task == "bottleneck" else [np.array([.5, 1.0, .1])]
    for _ in range(restarts - 1):
        starts.append(np.exp(rng.uniform(np.log([b[0] for b in bounds]), np.log([b[1] for b in bounds]))))

    def objective(log_parameters):
        expected = expected_spectrum(np.exp(log_parameters), sample_sizes, task)
        if np.any(~np.isfinite(expected.compressed())) or np.any(expected.compressed() <= 0):
            return 1e100
        likelihood = float(moments.Inference.ll_multinom(expected, data))
        return -likelihood if np.isfinite(likelihood) else 1e100

    records = []
    for method, options in [
        ("L-BFGS-B", {"maxiter": 300, "ftol": 1e-10, "gtol": 1e-5, "eps": 1e-4}),
        ("Powell", {"maxiter": 300, "ftol": 1e-10, "xtol": 1e-5}),
    ]:
        for start in starts:
            result = minimize(objective, np.log(start), method=method,
                              bounds=[tuple(np.log(b)) for b in bounds], options=options)
            pars = np.exp(result.x)
            expected = expected_spectrum(pars, sample_sizes, task)
            theta = float(moments.Inference.optimal_sfs_scaling(expected, data))
            boundary = any(abs(np.log(p / lo)) < 1e-3 or abs(np.log(p / hi)) < 1e-3 for p, (lo, hi) in zip(pars, bounds))
            records.append({"method": method, "start": start.tolist(), "parameters": pars.tolist(),
                            "negative_log_likelihood": float(result.fun), "theta": theta,
                            "success": bool(result.success), "message": str(result.message),
                            "at_boundary": boundary, "evaluations": int(result.nfev)})
    finite = [r for r in records if np.isfinite(r["negative_log_likelihood"]) and r["negative_log_likelihood"] < 1e90]
    if not finite:
        raise RuntimeError("no finite optimization result")
    # A better boundary/non-converged solution cannot be hidden by selecting a
    # worse interior result. Preserve the global best as a failed diagnostic.
    best = min(finite, key=lambda r: r["negative_log_likelihood"])
    eligible = [r for r in finite if r["success"] and not r["at_boundary"]]
    agreeing = [r for r in eligible if abs(r["negative_log_likelihood"] - best["negative_log_likelihood"]) < .01
                and np.max(np.abs(np.log(np.array(r["parameters"]) / best["parameters"]))) < .01]
    agreement = {"methods": sorted({r["method"] for r in agreeing}), "log_parameter_tolerance": .01,
                 "log_likelihood_tolerance": .01,
                 "passed": best["success"] and not best["at_boundary"] and {r["method"] for r in agreeing} == {"L-BFGS-B", "Powell"}}
    return {"task": task, **physical_parameters(task, best["parameters"], best["theta"], callable_length),
            "best_fit": best, "restarts": records, "agreeing_restarts": len(agreeing), "optimizer_agreement": agreement, "callable_length": callable_length,
            "segregating_sites_expected": float(data.sum()), "scientifically_ready": False,
            "release_status": "development_only_independent_full_validation_pending"}



def realsfs_diagnostics(log: str, tolerance: float = 1e-8) -> dict:
    changes = re.findall(r"lik\[([0-9]+)\]=[-+0-9.e]+ diff=([-+0-9.e]+)", log)
    norms = re.findall(r"Breaking EM\(s[rq]2\) at iter:([0-9]+), sqrt\(s[rq]2\):([-+0-9.e]+)", log)
    # ANGSD uses either likelihood change or the norm of a probability update.
    norm = float(norms[-1][1]) if norms else None
    change = float(changes[-1][1]) if changes else None
    return {"last_iteration": int(norms[-1][0]) if norms else (int(changes[-1][0]) if changes else None),
            "last_likelihood_change": change, "probability_update_norm": norm, "tolerance": tolerance,
            "convergence_observed": bool((change is not None and abs(change) < tolerance)
                                         or (norm is not None and abs(norm) < tolerance))}


def run_pipeline(directory: Path, reference: Path, task: str, output: Path, angsd: Path, realsfs: Path,
                 *, threads: int = 4, expected_length: int = SEQUENCE_LENGTH, restarts: int = 8, genotype_likelihood_model: int = 1) -> dict:
    """Execute the actual likelihood tools and retain every exact command/log."""
    if genotype_likelihood_model not in (1, 2):
        raise ValueError("likelihood model must be ANGSD1 or2")
    if threads < 1:
        raise ValueError("positive thread count required")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    audit = audit_inputs(directory, reference, task, expected_length)
    output.mkdir(parents=True)
    write_json(output / "inputs.json", audit)
    protocol = {"method": f"angsd_gl{genotype_likelihood_model}_saf1_realsfs_global_fold_moments_v5", "task": task,
                "source_sha256": sha256(Path(__file__)), "angsd_sha256": sha256(angsd),
                "realsfs_sha256": sha256(realsfs), "commands": [], "hard_genotype_calls": False,
                "snp_ascertainment": False, "ancestral_reference_assumed": False, "input_origin": "caller_supplied"}
    indices = []
    for population_index, population in enumerate(audit["population_order"]):
        listing = output / f"population_{population_index}.bams"
        listing.write_text("".join(s["path"] + "\n" for s in audit["samples"].values() if s["population"] == population))
        prefix = output / f"population_{population_index}"
        command = [str(angsd), "-b", str(listing), "-ref", str(reference.resolve()), "-anc", str(reference.resolve()),
                   "-GL", str(genotype_likelihood_model), "-doSaf", "1", "-minMapQ", "30", "-minQ", "20", "-P", str(threads), "-out", str(prefix)]
        protocol["commands"].append(command)
        write_json(output / "protocol.json", protocol)
        with (output / f"population_{population_index}.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        index = Path(str(prefix) + ".saf.idx")
        if not index.is_file() or index.stat().st_size == 0:
            raise RuntimeError("ANGSD produced no SAF index; inspect the population log")
        indices.append(str(index))
    command = [str(realsfs), *indices, "-fold", "1", "-P", str(threads), "-maxIter", "1000", "-tole", "1e-8", "-seed", "20261120"]
    protocol["commands"].append(command)
    write_json(output / "protocol.json", protocol)
    with (output / "folded.sfs").open("w") as result, (output / "realsfs.log").open("w") as log:
        subprocess.run(command, stdout=result, stderr=log, check=True)
    values = np.loadtxt(output / "folded.sfs")
    if values.ndim != 1 or not 0 < values.sum() <= expected_length + 1:
        raise ValueError("unexpected SFS dimensions or callable site total")
    estimate = fit_spectrum(values, audit["sample_sizes"], task, restarts=restarts)
    estimate["realsfs_diagnostics"] = realsfs_diagnostics((output / "realsfs.log").read_text())
    estimate["numerically_valid"] = estimate["optimizer_agreement"]["passed"] and estimate["realsfs_diagnostics"]["convergence_observed"]
    estimate["callable_fraction"] = float(values.sum() / expected_length)
    estimate["input_audit_sha256"] = sha256(output / "inputs.json")
    estimate["protocol_sha256"] = sha256(output / "protocol.json")
    estimate["sfs_sha256"] = sha256(output / "folded.sfs")
    write_json(output / "estimate.json", estimate)
    return estimate
