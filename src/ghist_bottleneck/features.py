from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from .arg_features import ARG_TIME_BINS


EPS = 1e-8
WINDOW_FEATURES = (
    "log1p_segregating_sites",
    "pi_per_site",
    "theta_w_per_site",
    "tajima_d",
    "singleton_fraction",
    "rare_fraction",
    "mean_folded_frequency",
    "heterozygosity",
)
LD_BINS = 12


@dataclass(frozen=True)
class FeatureBundle:
    windows: np.ndarray
    sfs: np.ndarray
    n_haplotypes: int
    sequence_length: float
    n_variants: int

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            windows=self.windows.astype(np.float32),
            sfs=self.sfs.astype(np.float32),
            n_haplotypes=np.int64(self.n_haplotypes),
            sequence_length=np.float64(self.sequence_length),
            n_variants=np.int64(self.n_variants),
        )

    @classmethod
    def load(cls, path: str | Path) -> "FeatureBundle":
        d = np.load(path)
        return cls(
            windows=d["windows"],
            sfs=d["sfs"],
            n_haplotypes=int(d["n_haplotypes"]),
            sequence_length=float(d["sequence_length"]),
            n_variants=int(d["n_variants"]),
        )


def _tajima_constants(n: int) -> tuple[float, float, float]:
    if n < 3:
        return 1.0, 0.0, 0.0
    a1 = np.sum(1.0 / np.arange(1, n))
    a2 = np.sum(1.0 / np.square(np.arange(1, n)))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return float(a1), float(e1), float(e2)


def summarize_variants(
    positions: np.ndarray,
    allele_counts: np.ndarray,
    called_haplotypes: np.ndarray,
    heterozygosity: np.ndarray,
    *,
    sequence_length: float,
    n_haplotypes: int,
    n_windows: int = 64,
    sfs_bins: int = 64,
    ld_positions: np.ndarray | None = None,
    ld_genotypes: np.ndarray | None = None,
    arg_hazards: np.ndarray | None = None,
) -> FeatureBundle:
    """Make fixed-shape, polarization-free summaries from variant arrays."""
    positions = np.asarray(positions, dtype=np.float64)
    ac = np.asarray(allele_counts, dtype=np.float64)
    called = np.asarray(called_haplotypes, dtype=np.float64)
    het = np.asarray(heterozygosity, dtype=np.float64)
    valid = (called >= 2) & (ac > 0) & (ac < called)
    positions, ac, called, het = positions[valid], ac[valid], called[valid], het[valid]
    folded = np.minimum(ac, called - ac)
    folded_freq = folded / called

    hist, _ = np.histogram(folded_freq, bins=np.linspace(0, 0.5, sfs_bins + 1))
    spectrum = hist.astype(np.float64)
    spectrum = np.log1p(spectrum)
    spectrum /= np.linalg.norm(spectrum) + EPS
    ld = _ld_decay(ld_positions, ld_genotypes, sequence_length)
    if arg_hazards is None:
        arg_context = np.zeros(ARG_TIME_BINS + 1, dtype=np.float64)
    else:
        hazards = np.asarray(arg_hazards, dtype=np.float64)
        if hazards.shape != (ARG_TIME_BINS,):
            raise ValueError(f"arg_hazards must have shape ({ARG_TIME_BINS},)")
        arg_context = np.concatenate((hazards, [1.0]))
    sfs = np.concatenate((spectrum, ld, arg_context))

    edges = np.linspace(0, sequence_length, n_windows + 1)
    win = np.zeros((n_windows, len(WINDOW_FEATURES)), dtype=np.float64)
    a1, e1, e2 = _tajima_constants(n_haplotypes)
    win_len = sequence_length / n_windows
    idx = np.clip(np.searchsorted(edges, positions, side="right") - 1, 0, n_windows - 1)
    for w in range(n_windows):
        m = idx == w
        s = int(m.sum())
        if not s:
            continue
        p = ac[m] / called[m]
        pi_total = np.sum(2 * p * (1 - p) * called[m] / np.maximum(called[m] - 1, 1))
        theta_total = s / max(a1, EPS)
        denom = np.sqrt(max(e1 * s + e2 * s * (s - 1), EPS))
        tajima = (pi_total - theta_total) / denom
        ff = folded_freq[m]
        win[w] = (
            # Express density per Mb so simulations can use shorter chromosomes
            # without creating a train/test scale cue.
            np.log1p(s * 1_000_000.0 / win_len),
            pi_total / win_len,
            theta_total / win_len,
            tajima,
            np.mean(folded[m] <= 1.01),
            np.mean(ff <= 0.1),
            np.mean(ff),
            np.mean(het[m]),
        )

    return FeatureBundle(win, sfs, n_haplotypes, sequence_length, len(positions))


def _ld_decay(
    positions: np.ndarray | None,
    genotypes: np.ndarray | None,
    sequence_length: float,
    *,
    n_bins: int = LD_BINS,
    pairs_per_bin: int = 512,
) -> np.ndarray:
    """Estimate dosage r² in log-spaced distance bins deterministically."""
    if positions is None or genotypes is None or len(positions) < 3:
        return np.zeros(n_bins, dtype=np.float64)
    positions = np.asarray(positions, dtype=np.float64)
    genotypes = np.asarray(genotypes, dtype=np.float64)
    order = np.argsort(positions)
    positions, genotypes = positions[order], genotypes[order]
    max_distance = max(200.0, min(sequence_length / 2, 5_000_000.0))
    edges = np.geomspace(50.0, max_distance, n_bins + 1)
    rng = np.random.default_rng(7081)
    result = np.zeros(n_bins, dtype=np.float64)
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        values = []
        for _ in range(pairs_per_bin * 3):
            i = int(rng.integers(0, len(positions) - 1))
            target = positions[i] + np.exp(rng.uniform(np.log(lo), np.log(hi)))
            j = int(np.searchsorted(positions, target))
            if j >= len(positions):
                continue
            distance = positions[j] - positions[i]
            if not (lo <= distance < hi):
                continue
            a, b = genotypes[i], genotypes[j]
            valid = np.isfinite(a) & np.isfinite(b)
            if valid.sum() < 4:
                continue
            aa, bb = a[valid], b[valid]
            aa, bb = aa - aa.mean(), bb - bb.mean()
            denom = np.sqrt(np.sum(aa * aa) * np.sum(bb * bb))
            if denom > EPS:
                values.append(float(np.square(np.sum(aa * bb) / denom)))
            if len(values) >= pairs_per_bin:
                break
        if values:
            result[k] = np.mean(values)
    return result


def extract_vcf(path: str | Path, *, n_windows: int = 64, sfs_bins: int = 64) -> FeatureBundle:
    vcf = pysam.VariantFile(str(path))
    samples = list(vcf.header.samples)
    positions: list[float] = []
    allele_counts: list[int] = []
    called_haplotypes: list[int] = []
    heterozygosity: list[float] = []
    max_pos = 0
    inferred_ploidy = 2
    # Bounded deterministic reservoir for LD; avoids retaining the full matrix.
    ld_capacity = 5000
    ld_positions: list[float] = []
    ld_genotypes: list[np.ndarray] = []
    rng = np.random.default_rng(7081)
    eligible_seen = 0
    for rec in vcf:
        max_pos = max(max_pos, rec.pos)
        if len(rec.alleles or ()) != 2 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        ac = called = het = called_genotypes = 0
        dosages = np.full(len(samples), np.nan, dtype=np.float32)
        for sample_i, sample in enumerate(samples):
            gt = rec.samples[sample].get("GT")
            if gt is None:
                continue
            vals = [a for a in gt if a is not None and a >= 0]
            if not vals:
                continue
            inferred_ploidy = max(inferred_ploidy, len(vals))
            called += len(vals)
            ac += sum(a == 1 for a in vals)
            dosages[sample_i] = sum(a == 1 for a in vals)
            if len(vals) == 2:
                het += vals[0] != vals[1]
                called_genotypes += 1
        if 0 < ac < called:
            positions.append(rec.pos - 1)
            allele_counts.append(ac)
            called_haplotypes.append(called)
            heterozygosity.append(het / max(called_genotypes, 1))
            eligible_seen += 1
            if len(ld_positions) < ld_capacity:
                ld_positions.append(rec.pos - 1)
                ld_genotypes.append(dosages)
            else:
                replace = int(rng.integers(0, eligible_seen))
                if replace < ld_capacity:
                    ld_positions[replace] = rec.pos - 1
                    ld_genotypes[replace] = dosages
    seq_len = max(float(max_pos), 100_000_000.0)
    return summarize_variants(
        np.asarray(positions),
        np.asarray(allele_counts),
        np.asarray(called_haplotypes),
        np.asarray(heterozygosity),
        sequence_length=seq_len,
        n_haplotypes=len(samples) * inferred_ploidy,
        n_windows=n_windows,
        sfs_bins=sfs_bins,
        ld_positions=np.asarray(ld_positions),
        ld_genotypes=np.asarray(ld_genotypes),
    )
