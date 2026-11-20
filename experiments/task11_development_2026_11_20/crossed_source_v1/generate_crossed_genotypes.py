"""Independent/crossed nuisance assignments; original genotype mechanics retained."""
import gzip
from pathlib import Path

import numpy as np
from scripts.develop_task11_pedigree import (
    LABELS,
    dump,
    expand_inbreeding,
    focal_pedigree,
    gamete,
    read_map,
)


def validate_assignment(config):
    if len(config["ancestry_kinds"]) != len(LABELS) or not set(config["ancestry_kinds"]) <= {"unrelated", "first_cousins", "siblings"}:
        raise ValueError("Explicit ancestry kind required for every focal family")
    if len(config["third_branches"]) != len(LABELS) or any(type(value) is not bool for value in config["third_branches"]):
        raise ValueError("Explicit third-branch booleans required for every focal family")
    if not 1 <= config["min_branch_size"] <= config["max_branch_size"] <= 4:
        raise ValueError("Unsupported branch size range")
    if config["complex_topology"] and any(not config["third_branches"][i] for i, label in enumerate(LABELS) if label.startswith("avuncular")):
        raise ValueError("Mixed avuncular topology requires its third branch")


def generate_genotypes(map_path: Path, out: Path, config: dict) -> dict:
    validate_assignment(config)
    import msprime

    rng = np.random.default_rng(config["seed"])
    maps = read_map(map_path)
    records = []
    samples, families, truth_pairs = [], [], []
    next_sample = 0
    for index, label in enumerate(LABELS):
        sizes = tuple(int(x) for x in rng.integers(config["min_branch_size"], config["max_branch_size"] + 1, 2))
        complex_topology = bool(config.get("complex_topology", False))
        third = config["third_branches"][index]
        pedigree = focal_pedigree(label, sizes, third, complex_topology)
        parents = expand_inbreeding(pedigree, config["ancestry_kinds"][index])
        names = {n: f"synthetic11_{next_sample + i}" for i, n in enumerate(pedigree["observed"])}
        next_sample += len(names)
        families.append(list(names.values()))
        samples.extend(names.values())
        truth_pairs.extend({"indID1": names[a], "indID2": names[b], "relation": relation,
                            "family": index} for a, b, relation in pedigree["pairs"])
        records.append({"parents": parents, "observed": pedigree["observed"],
                        "kind": config["ancestry_kinds"][index], "names": names})
    founder_count = sum(sum(pair is None for pair in r["parents"].values()) for r in records)
    hbd_cm = np.zeros(len(samples))
    genome_cm, markers, missing_calls = 0.0, 0, 0
    with gzip.open(out / "synthetic.vcf.gz", "wt") as destination:
        destination.write("##fileformat=VCFv4.2\n")
        destination.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Known simulated phase">\n')
        for chrom, (bp, _) in maps.items():
            destination.write(f"##contig=<ID={chrom},length={int(bp[-1]) + 1}>\n")
        destination.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for chrom, (map_bp, map_cm) in maps.items():
            positions = np.concatenate(([0.0], map_bp.astype(float), [float(map_bp[-1] + 1)]))
            # Remove duplicate 0 if the original map begins there.
            positions = np.unique(positions)
            cm_boundaries = np.interp(positions, map_bp, map_cm)
            rates = np.diff(cm_boundaries) / 100.0 / np.diff(positions)
            ts = msprime.sim_ancestry(samples=founder_count, population_size=config["population_size"],
                recombination_rate=msprime.RateMap(position=positions, rate=rates),
                random_seed=int(rng.integers(1, 2**31)))
            ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=int(rng.integers(1, 2**31)))
            site_bp, alleles = [], []
            last = -500
            for variant in ts.variants():
                bp = int(variant.site.position)
                g = variant.genotypes
                if (bp < map_bp[0] or bp > map_bp[-1] or bp - last < 500
                    or len(variant.alleles) != 2 or min(int(g.sum()), len(g) - int(g.sum())) < 2):
                    continue
                site_bp.append(bp)
                alleles.append(g.copy())
                last = bp
            if len(site_bp) < 100:
                raise ValueError(f"Insufficient segregating markers on chromosome {chrom}")
            bp = np.array(site_bp)
            cm = np.interp(bp, map_bp, map_cm)
            founder_alleles = np.asarray(alleles, dtype=np.int8).T
            offset, ancestry_samples = 0, []
            for record in records:
                ancestry = {}
                for name, pair in record["parents"].items():
                    if pair is None:
                        ancestry[name] = np.repeat(np.arange(offset, offset + 2)[:, None], len(bp), axis=1)
                        offset += 2
                    else:
                        ancestry[name] = np.stack([gamete(ancestry[p], cm, rng) for p in pair])
                ancestry_samples.extend(ancestry[n] for n in record["observed"])
            ancestors = np.stack(ancestry_samples)
            genotypes = founder_alleles[ancestors, np.arange(len(bp))]
            is_hbd = ancestors[:, 0] == ancestors[:, 1]
            hbd_cm += (is_hbd[:, :-1] * np.diff(cm)).sum(axis=1)
            genome_cm += float(map_cm[-1] - map_cm[0])
            missing = rng.random((len(samples), len(bp))) < config["missing_rate"]
            missing_calls += int(missing.sum())
            for j, position in enumerate(bp):
                calls = [".|." if missing[i, j] else f"{int(g[0])}|{int(g[1])}"
                         for i, g in enumerate(genotypes[:, :, j])]
                destination.write(f"{chrom}\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
            markers += len(bp)
    truth = {"samples": samples, "families": families, "pairs": truth_pairs,
             "pedigrees": records, "genome_cm": genome_cm, "markers": markers,
             "founder_individuals": founder_count, "seed": config["seed"],
             "missing_genotypes": missing_calls, "true_hbd_cm": dict(zip(samples, hbd_cm.tolist()))}
    dump(out / "truth.json", truth)
    return truth
