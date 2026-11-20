# Readiness follow-up, 20 November 2026

Task 9 now has an explicitly selected, versioned release candidate that passes
its scientific and production checks in the tested exact-IBD input domain.
This is the strongest working challenge pipeline. Tasks 1–8 and 10–12 still
have scientific validation or accuracy gaps. The two later low-coverage
challenges remain unimplemented. No current final-data inference or competition
submission has been performed in this follow-up.

The user authorized direct merging and parallel scientific review, while asking
that final-data runs wait. Three reviewers cover demography, sweeps and
relatedness; the coordinating reviewer covers Task 1 and integration. This
report supplements the earlier [scientific review](scientific_validation_2026_11_20.md),
whose archived results and historical artifacts remain preserved.

## What is working

| Challenge | Current evidence | Status |
|---|---|---|
| 1 — bottleneck | Software and input-contract repairs pass; new 24-case, 100 Mb panel running | Scientific confirmation pending |
| 2 — secondary contact | All original cases now finish; hidden-history and contact-time errors remain | Runs, not scientifically ready |
| 3 — growth/BGS | New 10 Mb neutral simulations recover strong growth; weak growth remains poorly identified | Runs, broader validation needed |
| 4 — split/migration/BGS | Numerical stall fixed; Q1 fit finishes, but split-time error is +47% | Runs, accuracy/calibration gaps |
| 5 — hard sweep | Production VCF-to-BED path implemented; fresh full-chromosome panel running | Scientific confirmation pending |
| 6 — soft sweep | Eight new development configurations fail recovery criteria | Not ready |
| 7 — sweep/BGS | Input handling repaired; exact-estimator BGS validation missing | Not ready |
| 8 — multiple sweeps/BGS | Existing detection/localization failures remain | Not ready |
| 9 — exact IBD relationships | Fixed candidate passes scientific and production gates | Validated candidate in tested domain |
| 10 — human pedigrees | Separate complete-pipeline review underway | Not yet certified |
| 11 — dog pedigrees | Separate complete-pipeline review found consensus accuracy failures | Not ready |
| 12 — background-related families | Corrected independent-reference candidate fails improvement gate | Not ready |

**Task 9 — exact IBD relationships.** The ten-label candidate passed eight
executions of the real predictor and independent submission validator on two
untouched simulated datasets containing 20,023 unique pairs. Exact accuracy is
89.10% under the matched crossover model and 87.20% under a Poisson crossover
shift. Reversing pairs, reordering rows/columns and adding a UTF-8 BOM leaves
predictions unchanged. Model, map, feature contract and runtime checksums guard
release selection. These are conditional accuracy estimates, not a claim that
every relationship is correct. See the [Task 9 report](task9_readiness_2026_11_20.md).

The official Task 9 input consists of exact IBD segments. Genotype simulation,
phasing and IBD detector validation are therefore not prerequisites for this
task; those remain essential for the complete Tasks 10–12 pipelines. Pair
coverage is checked against all unordered pairs represented in the segment
CSV. The absence of an independent zero-IBD pair roster is documented explicitly.

## Repairs and continuing validation

**Tasks 2 and 4 — bounded numerical integration.** The installed moments solver
can repeatedly restore an invalid spectrum at its minimum time step without
advancing simulation time. Instrumented reproductions establish this for both
tasks. A bounded sequence of complete non-adaptive integrations now checks
valid likelihood bins and agreement between time resolutions. It does not clip
negative entries, discard observed bins, or report a fallback answer. Optimizer
termination is recorded separately from finite likelihood and restart agreement.
Atomic spectrum and fit checkpoints preserve expensive simulation work and
validate configuration and content before reuse. Newly allocated random seed
streams cannot collide; legacy allocation is restricted to reproducing the
original audit panel.

The previously stalled unscaled Task 4 fit completed in 1 minute 52 seconds
after this repair. The same preserved tree and physical model were used.
Remaining errors include roughly +47% split time, −16% EUR population size and
−23% migration. Execution correctness does not establish model adequacy.
All 22 original neutral cases and 28 new 10 Mb cases completed. Eleven selected
fits have abnormal optimizer termination; release guards reject these fits.
Strong-growth recovery improves at the longer length, but weak exponential
growth still has up to 71.6% timing error. See the
[demographic report](readiness_demography_2026_11_20.md) for every denominator,
failure stratum and numerical sensitivity check.

**Task 1 — full-length recovery and input checks.** A prespecified 24-case panel
uses 100 Mb chromosomes, 57 diploid individuals, six diagnostic parameter
points, two independent replicates and matched versus perturbed ancestral
histories. All cases use random optimizer initialization. Independent dadi
diffusion fits accompany the matched moments fits. Per-case block spectra are
saved before fitting; failures and incomplete cases remain visible. Two cases
also exercise the actual VCF and block-bootstrap readers, requiring exact
agreement with tree-sequence spectra. This small diagnostic panel does not
establish full-prior accuracy or confidence-interval coverage.

Spectrum projection now groups repeated allele counts and directly counts sites
at the requested sample size. This preserves exact hypergeometric projection
while avoiding one distribution evaluation per variant. Unpolarized VCF input
requires folding rather than treating ALT as the ancestral-derived orientation. Fits expose actual
optimizer termination. Bootstrap fits accept valid boundary reference values,
reject mismatched mutation rates and exposure, retain failed replicate records,
verify the reference VCF and spectrum identity, and omit intervals if any
replicate fails or its selected optimizer does not converge. Percentile intervals are explicitly
marked as having unvalidated coverage.

**Tasks 5–8 — consistent accessibility and blind localization.** The production
hard-sweep path applies accessibility to both scanning and fine localization.
Inaccessible gaps cannot act as sweep evidence; normalization excludes unusable
windows. BED containment follows half-open coordinates. A prespecified fresh
whole-chromosome Task 5 protocol checks production VCF reading, localization,
BED writing and the submission validator. Its resource pilot is separate from
accuracy-based selection. The old Task 6 diagnostic set is explicitly development
data for evaluating SNP-window H12 methods; it cannot be reused as fresh
confirmation. Only Task 8 explicitly specifies bottleneck/growth in addition to
background selection; the Task 7 gate now states its matched-BGS requirement
without attributing unrequired demographic history to the specification.

**Task 12 — independent reference validation.** A fixed corrected refit was
tested on four new pedigree seeds, using calibration-only correction references.
All 40 CSVs passed pair coverage and format checks, but exact accuracy fell from
37.64% to 37.37%. All four seed-average changes were negative. The candidate is
not promoted; the component experiment does not certify joint family
reconstruction. See the [Task 12 report](task12_readiness_2026_11_20.md).

## Scientific limits

Task 5's earlier 25.4 kb median localization error was retrospective and is not
fresh release evidence. Tasks 6 and 8 fail important localization/detection
checks. Tasks 1–4 require broader model recovery, artifact regeneration where
provenance was invalidated, and uncertainty validation. Tasks 7 and 10–12 lack
complete matched pipeline validation; pair-template performance is insufficient
to certify family reconstruction.

During this follow-up, two reviewers inadvertently encountered embedded public
testing-answer sections while retrieving official specifications. They reported
the exposure immediately. Answer values are excluded from the new evidence and
were not used to choose targets, parameters, candidate methods or thresholds.
The Task 9 candidate and its two confirmation simulations were already fixed
before this exposure. New sweep designs use mechanical random positions and
prespecified cells, but this follow-up cannot be described as wholly blinded
to public answers. The earlier frozen scientific evidence is unchanged.

Claude's requested duplicate reviews remain incomplete because its account
quota rejects execution. Earlier interrupted attempts do not count as completed
reviews. This limitation is separate from the completed parallel reviews and
measured experimental evidence.

The isolated source snapshot passes all **340 tests** on macOS/Python 3.12,
Linux ARM/Python 3.12 and Betty/Linux x86-64/Python 3.11. Ruff also passes.
The snapshot excludes concurrent work from the separate submission and
pedigree-pipeline chats, whose integration is tracked separately. Cross-review
reproductions added guards for incomplete or duplicate projection sets,
nonfinite freeze values, unknown checkpoint status, unverified metric records,
invalid sweep markers and mismatched download-manifest field names. Initial
remote harness failures and their corrections are preserved in the verification
records rather than omitted.

Compact protocols, checkpoints, job scripts and verification records are under
[`experiments/readiness_2026_11_20`](../experiments/readiness_2026_11_20/).
Scientific accuracy claims apply only to the stated experiments. No finite test
suite or simulation panel can establish that software is bug-free or that every
inference will be correct.
