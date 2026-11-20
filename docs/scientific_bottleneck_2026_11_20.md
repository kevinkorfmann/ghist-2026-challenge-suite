# Bottleneck scientific validation — November 20, 2026

The direct two-epoch spectrum estimator runs, but this review does not establish
uniform accuracy or interval coverage. No public challenge answers were used.
No final data were fitted and no submission was made.

## Scientific corrections

The three bottleneck simulation entry points used global `ploidy=1` to request
haplotype counts while treating ancestral population size as diploid effective
size. In msprime, global ploidy also scales the population's coalescent time.
The old nominal size prior therefore corresponded to half the intended diploid
size. Event times themselves were still in generations; this is not a blanket
factor-of-two correction to every fitted time. See the
[msprime ploidy documentation](https://tskit.dev/msprime/docs/stable/ancestry.html#ploidy).

Simulation now uses gene-copy SampleSets with global diploid ploidy. New datasets,
merged shards, and learned checkpoints record `coalescent_ploidy=2`. Training,
learned inference, chromosome validation and multiscale calibration reject legacy
or unverified scaling. Old artifacts are preserved for audit and require
regeneration/retraining to certify the intended prior. This does not change the
direct moments fit to an observed spectrum.

Training also refuses nonempty output directories, and inference checks every
network member, reference estimator and chromosome calibration separately. A
new normalization file therefore cannot make leftover legacy weights eligible.

The independent dadi reference failed on every new nominal case because its
wrapper passed `iprint` to SciPy 1.18. The reference now calls the same public
log-space L-BFGS-B optimizer directly, preserving the bounds, finite-difference
step and function-evaluation budget. The current
[SciPy API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html)
does not accept `iprint`. A deterministic known-spectrum recovery test exercises
the real diffusion fit. The pre-fix failure record is retained alongside the
successful refit of the exact same simulated spectra.

Both reference methods now reject expected spectra that would silently exclude
observed likelihood bins. Variant-only VCF fitting with incomplete calls requires
independently measured callable exposure for the projection, exposed as
`--callable-sequence-length`. A variant-only VCF cannot establish the number of
callable invariant sites.

## Independent recovery panel

[Protocol and results](../experiments/scientific_2026_11_20/bottleneck/recovery.json)
contain 24 independent 10 Mb chromosomes, 114 gene copies, two replicates at each
of six prespecified parameter points, and two histories. The matched history is
a single decline. The misspecified history also has an ancestral size epoch
from four to eight times the decline time. The fit never receives the simulated
truth as an initial estimate. No estimator is selected or calibrated on this
panel. Rate, seed, package and source provenance are recorded.
The recorded two-replicate panel has unique seeds. A subsequent cross-review
found that the original arithmetic seed allocation collided at eleven or more
replicates; new panels use disjoint seed triples, and reused panels are checked
for complete case identities and seed uniqueness.

| History and method | Successful fits | Median absolute relative timing error | Median absolute relative fraction error |
|---|---:|---:|---:|
| Matched, moments | 12/12 | 19.03% | 1.71% |
| Matched, dadi | 12/12 | 19.21% | 1.88% |
| Extra ancestral epoch, moments | 12/12 | 35.90% | 14.86% |

These medians conceal severe failures for weak declines. At a true fraction of
0.94 and time 4,000 generations, one replicate produces approximately 163,973
generations with moments and 148,782 with dadi. The largest matched timing
relative error is 39.99 (about 3,999%), despite finite optimizer outputs.
Agreement between diffusion and moment solvers therefore cannot establish
biological accuracy: they share the same demographic model and data.
The [descriptive summary](../experiments/scientific_2026_11_20/bottleneck/summary.json)
also reports likelihood gains over a constant-size model; these are not
calibrated hypothesis tests or confidence intervals.

The panel is a diagnostic stress test, not a draw from the full challenge prior
and not a target-length 100 Mb benchmark. Two replicates per point cannot
establish coverage or a reliable tail-failure frequency. Historical chromosome
validation is informative about its actual, shifted size distribution, but
cannot certify the originally claimed diploid prior.

## Remaining evidence gates

Regenerate the intended-prior training/calibration and chromosome panels;
evaluate target-length recovery and block-interval coverage with enough
independent chromosomes; assess linkage-based references on the same held-out
datasets; and quantify missingness, rate errors and hidden-history shifts.
At a fixed SFS, inferred size and time scale inversely with assumed mutation
rate, so numerical convergence does not remove mutation-rate uncertainty.
Final-input provenance and a fresh final-data run remain deferred.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.validate_bottleneck_science \
  --spectra-from experiments/scientific_2026_11_20/bottleneck/recovery.json \
  --out artifacts/bottleneck_scientific_refit.json
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.summarize_bottleneck_science \
  --input experiments/scientific_2026_11_20/bottleneck/recovery.json \
  --out experiments/scientific_2026_11_20/bottleneck/summary.json
```

Omitting `--spectra-from` generates a new panel with the corrected disjoint seed
allocation, rather than recreating the archived arithmetic seed sequence.
