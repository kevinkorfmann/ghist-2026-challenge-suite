# Independent scientific validation: Tasks 2–4

This review found scientific implementation errors that invalidate parts of the
historical calibration evidence. Those errors are corrected, and old artifacts
cannot silently pass the affected current inference or freeze entry points.
None of Tasks 2–4 is certified ready for final submission. No public answers were
read for this review, no parameter was tuned to challenge answers, and no final
data or submission service was used.

## Corrections and their scientific consequences

**Task 2 fixed differences.** Training tree-sequence extraction and two
calibration simulators discarded all four joint-spectrum corners. Only the two
globally monomorphic corners should be discarded: `(0, n2)` and `(n1, 0)` are
fixed differences between populations, informative about divergence and gene
flow. The observed VCF scanner already retained these sites. Thus the defect
created an observation/simulation mismatch in addition to losing information.
The extraction now retains fixed differences, with a deterministic tskit
regression covering both orientations and a globally fixed site. Existing
learned checkpoints and bootstrap/calibration evidence must be regenerated.

The hidden-history simulator also omitted a mainland size change whenever it
fell after the onset of contact; the corresponding diffusion model retained it.
The simulator now follows the stated history at all valid contact fractions.
The historical hidden-epoch convention is still an assumption, not evidence
about the unknown current challenge generator.

Analytical spectrum-surrogate training previously dropped the model mask and
clipped negative numerical results to zero. It now excludes masked monomorphic
bins, retries unstable integration adaptively, and rejects an invalid spectrum.
The trained spectrum transform is recorded and used by tabular prediction;
training rejects mismatched transforms. A version marker propagates through
datasets, models, calibration summaries, and freeze checks, so old evidence
cannot be relabeled automatically as corrected validation.

**Tasks 3 and 4 SLiM scaling.** Previous recovery drivers multiplied population
sizes and times by computational scaling factor Q, divided mutation rate by Q
(and migration by Q for Task 4), and then asked stdpopsim to apply Q scaling.
Recombination maps and selection coefficients were left in physical units.
Consequently the model preserved `N*mu` and demographic time ratios, but targeted
Q-fold larger `N*r` and `N*s` than the stated analysis model. Good recovery in
those simulations did not establish matched background-selection recovery.

The corrected drivers pass physical sizes, times, mutation, migration, map, and
DFE parameters directly to stdpopsim; the engine alone applies computational
rescaling. This follows the [stdpopsim scaling API](https://popsim-consortium.github.io/stdpopsim-docs/stable/api.html)
and [scaling tutorial](https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html).
All relevant aggregation, calibration, and freeze scripts reject legacy scaling
records. Task 4 recovery fits no longer receive the generating truth as optimizer
initialization, and its aggregate fitter uses the recorded mutation rate.

**Likelihood domain in Tasks 2–4.** A finite moments likelihood can omit bins
whose numerical predictions are negative, nonfinite, zero despite positive
observations, or additionally masked. The installed moments `Inference.py`
explicitly uses masked logarithms and warnings for these cases. Model callbacks
now keep the observed likelihood domain fixed, retry stiff two-population
integrations adaptively, and reject invalid trajectories. Independent tests
exercise each failure mode. Chromosome growth comparison is explicitly labeled
as conditional spectrum-shape likelihood, because its held-out amplitude is
profiled rather than predicted.

## Independent evidence

The reproducible records are in
[`experiments/scientific_2026_11_20/demography`](../experiments/scientific_2026_11_20/demography/).
The neutral panel uses `scripts/validate_demography_independently.py`: two
replicate seeds per independently selected configuration, 2 Mb linked
chromosomes, msprime ancestry with diploid population-size conventions, and
infinite-sites mutations. Simulation truth is never passed to an optimizer.
Tasks 2 and 4 use the same generic scaled initialization at every truth point;
Task 3 uses global random restarts. These diagnostics with fixed initialization
do not establish global optimizer reliability. The initial unguarded run was
interrupted after the likelihood-domain defect was found, before any completed
result was used; the guarded panel uses unchanged simulation seeds and truths.

The panel contrasts simple and historical hidden-epoch secondary-contact
models, independent hidden-size changes and asymmetric migration, exponential
and sudden growth, and two split/migration projection sizes. Asymmetric cases
compare the symmetric estimator to the generating geometric-mean migration
rate; under misspecification there is no requirement that this pseudo-target be
identifiable. The generated [summary](../experiments/scientific_2026_11_20/demography/summary.json)
records every completed replicate and parameter error, including failures and
numerical diagnostics. Neutral runs are not BGS validation or whole-genome
uncertainty coverage.

The neutral panel was bounded at approximately 20 minutes of wall time on
eight CPU workers. Thirteen of 22 planned records completed: one Task 2, all
eight Task 3, and four Task 4 records. Seven Task 2 records and two Task 4
records remained incomplete when the run was stopped. They are explicitly
listed in the [execution record](../experiments/scientific_2026_11_20/demography/execution.json),
with their original configurations and seeds. The CPU-active fits
did not finish within this diagnostic budget; solver profiling is still needed.
Completed-only
results are a censored subset and must not be reported as an overall recovery
rate. The remaining records require bounded solver profiling and reruns before
they can contribute scientific evidence.

The panel already exposes limitations independently of BGS. At the weak-growth
point (1.6×, 3,100 generations), one matching exponential-model replicate has
approximately +355% timing error and −31% growth-ratio error. Under asymmetric
secondary-contact migration, a symmetric simple fit has approximately +158%
contact-time error in one replicate. These are retained outcomes, not discarded
outliers or targets for tuning.

Task 4 also shows replicate sensitivity without selection: the second matched
replicate at the first truth point underestimates split time by 37% and
migration by 42% at projection 20×16. Projection 12×10 changes those errors
to 37% and 46%, with a 34% overestimate of the smaller population size.
All completed fits pass the numerical optimizer-support diagnostic; that
diagnostic therefore cannot serve as a scientific accuracy certificate.

A separate mechanistic consistency experiment uses 10,000 independent unlinked
msprime genealogies for each of two early/late-contact histories. The corrected
historical diffusion model and empirical branch spectra agree to within 0.21%
in total branch length; normalized spectrum L1 differences are 1.94–2.05%.
Fixed-difference mass is present in both. This verifies implementation
consistency for the specified history, not the appropriateness of that history
for 2026.

The corrected BGS simulations use stdpopsim 0.3.0, SLiM 4.3, the Drosophila
Comeron map (halved for Task 4), exonic Gamma_H17 selection, chromosome 2L
5–6 Mb, and a 2N burn-in. Task 3 uses ancestral size 7,000, duration 9,200,
growth ratio 3.7, and 20 diploid samples. Task 4 uses ancestral size 11,000,
population sizes 8,000/19,000, split time 7,800, migration rate 0.00007,
and 12/10 diploid samples. These points were chosen independently of observed
challenge estimates.

| BGS diagnostic | Signed relative errors |
|---|---|
| Growth Q=1 | time +27.52%; ratio +4.59% |
| Growth Q=5 | time +8.29%; ratio −0.72% |
| Growth Q=10 | time −3.92%; ratio +6.34% |
| Split migration Q=5 | AFR −5.06%; EUR −10.48%; time −6.36%; migration −1.39% |
| Split migration Q=10 | AFR −27.27%; EUR −13.42%; time +16.05%; migration −5.89% |

The four Q=5/10 jobs completed as Slurm array **10849222**. Their recorded final
likelihoods were independently recomputed with the new likelihood-domain guard;
all included every observed bin and agreed within 1e-12 log units. Good
optimizer agreement did not prevent substantial biological estimation error.
On the exponential growth generator, Q=5 slightly favored exponential growth
(0.49 log units), while Q=10 favored sudden growth (2.09 log units). This is
evidence of weak discrimination in these small datasets, not grounds to choose
a convention. Cross-Q runs use independent seeds, so differences cannot be
attributed solely to the scaling approximation. A separate unscaled Q=1 array
is recorded as **10849487**. Its growth run also weakly favors the incorrect
sudden convention (0.95 log units), while timing recovery is appreciably worse
than in the Q=5/10 runs. Its results belong in the same evidence directory.

The Q=1 split-migration job produced its tree sequence but did not finish the
inference stage within a 15-minute diagnostic budget. It was cancelled at
900 seconds; no fit or recovery claim is reported for that run. The execution
record preserves its job status, source command, remote tree-sequence location,
and SHA-256 digest for later analysis. Its incompleteness also prevents using
the completed Q=5/10 fits to establish adequate scaling accuracy.

## Remaining gates

Task 2 needs regenerated simulation features, checkpoints, calibration and
interval evidence, independent recovery across current plausible hidden
histories, and likelihood/model discrepancy checks on current inputs. The
historical convention cannot be certified merely by matching its own simulator.

Task 3 needs replicated model discrimination and coverage under physical-unit
BGS, weak-growth and timing boundaries, long burn-in, selection/map shifts,
and full chromosome exposure. The existing chromosome holdouts do not
consistently identify one growth convention, and corrected small BGS runs do
not resolve that ambiguity.

Task 4 needs a regenerated broad physical-unit BGS grid, grouped calibration
selection with a separate untouched evaluation set, interval coverage, and
projection/chromosome stability under realistic selection. The present smokes
exceed several historical recovery tolerances. Two replicates or one region
cannot estimate reliable uncertainty coverage for any task.

Passing regression tests establishes that the corrected contracts and numerical
guards work. It does not establish accurate inference for unseen demographic
histories, current final data, or all parameter regimes.
