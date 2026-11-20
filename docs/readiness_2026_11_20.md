# Readiness follow-up, 20 November 2026

Task 9 now has an explicitly selected, versioned release candidate that passes
its scientific and production checks in the tested exact-IBD input domain.
This is the strongest working challenge pipeline. Tasks 1–8 and 10–12 still
have scientific validation or accuracy gaps. The two later low-coverage
challenges are being developed in a separate, unmerged continuation. No current
final-data inference or competition submission has been performed in this follow-up.

The user authorized direct merging and parallel scientific review, while asking
that final-data runs wait. Three reviewers cover demography, sweeps and
relatedness; the coordinating reviewer covers Task 1 and integration. This
report supplements the earlier [scientific review](scientific_validation_2026_11_20.md),
whose archived results and historical artifacts remain preserved.

## What is working

| Challenge | Current evidence | Status |
|---|---|---|
| 1 — bottleneck | Export checks pass; optimizer repair in development; completed weak-decline case has major recovery errors | Not ready; 24-case panel incomplete |
| 2 — secondary contact | All original cases now finish; hidden-history and contact-time errors remain | Runs, not scientifically ready |
| 3 — growth/BGS | New 10 Mb neutral simulations recover strong growth; weak growth remains poorly identified | Runs, broader validation needed |
| 4 — split/migration/BGS | Numerical stall fixed; Q1 fit finishes, but split-time error is +47% | Runs, accuracy/calibration gaps |
| 5 — hard sweep | Three completed cases: 2/3 interval coverage, 0/3 agreement passes; same-seed 22 GiB retry also exhausted memory | Not ready; simulation feasibility and recovery gaps |
| 6 — soft sweep | Eight new development configurations fail recovery criteria | Not ready |
| 7 — sweep/BGS | New masked full-input path passes real RAiSD software smoke; matching biological validation missing | Not ready |
| 8 — multiple sweeps/BGS | Existing detection/localization failures remain | Not ready |
| 9 — exact IBD relationships | Fixed candidate passes scientific and production gates | Validated candidate in tested domain |
| 10 — human pedigrees | Complete pipeline produces all 153 pairs in one family; only 58 exact labels correct | Runs, accuracy inadequate |
| 11 — dog pedigrees | Two cohorts: classifier 96/96 sibling pairs correct, final consensus 16/96 | Consensus fails; broader relationships unvalidated |
| 12 — background-related families | Full-histogram repair passes software checks but fails fresh accuracy criteria | Not ready |
| 13–14 — low-coverage BAM challenges | New implementation and synthetic validation in separate continuation | Development, not validated |

**Task 9 — exact IBD relationships.** The ten-label candidate passed eight
executions of the real predictor and independent submission validator on two
untouched simulated datasets containing 20,023 unique pairs. Exact accuracy is
89.10% under the matched crossover model and 87.20% under a Poisson crossover
shift. Reversing pairs, reordering rows/columns and adding a UTF-8 BOM leaves
predictions unchanged. Model, map, feature contract and runtime checksums guard
release selection. These are conditional accuracy estimates, not a claim that
every relationship is correct. See the [Task 9 report](task9_readiness_2026_11_20.md).

Four further production executions with the updated shared validator and
download-manifest-compatible predictor reproduce every prediction CSV byte for
byte. They reuse the same confirmation data and are compatibility evidence,
not additional independent biological replicates. The offline
[submission handoff](task9_submission_handoff_2026_11_20.md) keeps current final
input and prediction fields explicitly pending.

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
also exercised the actual VCF and block-bootstrap readers, with exact
agreement with tree-sequence spectra. The corrected export passed for both
100 Mb cases, including all 50 blocks. The old and new case 0 block spectra are
identical; the case 12 check binds the preserved tree to its original spectrum.
The [export proof](../experiments/readiness_2026_11_20/bottleneck/vcf_revalidation_v2/verification.json)
retains the exact validator source snapshot, separate from the original panel.
This small diagnostic panel does not establish full-prior accuracy or confidence-interval coverage.

The panel remains incomplete. Its completed weak-decline case 16 already shows
large physical-parameter errors in both solvers despite reported convergence.
Moments also flags poor restart agreement. This failure remains in the original
results; numerical optimizer repair cannot by itself establish identifiability
or reliable biological recovery.

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

Before inspecting Task 5 accuracy, resource measurements restricted this pass
to the first six balanced pilot cases. All 48 originally planned cases remain
in the accounting; the other 42 are explicitly not run. Three constant-size
cases completed with localization errors of 31.1, 39.9 and 108.0 kb, two covered
targets and no agreement-gate passes. The three initial variable-history
attempts failed operationally; one has a confirmed out-of-memory kill. A
bounded same-seed retry also exhausted its 22 GiB limit after 29 minutes 44 seconds;
Slurm and stderr independently confirm the OOM event. No further simulation
resource escalation is running.
The tiny completed subset cannot establish full-panel recovery. See the
[sweep report](readiness_sweeps_2026_11_20.md) for failures and resource policy.

**Tasks 10 and 11 — complete genotype-to-pedigree diagnostics.** A new
18-person human family passed through founder-genotype simulation, inheritance,
IBD detection, classification, reconstruction and consensus. All 153 pair
outputs were produced, but only 58 exact relationships and 117 meioses counts
were correct. One family does not supply 153 independent biological replicates.
For Task 11, two cohorts contain 16 families and 96 sibling pairs; clustering
and the flat sibling calls succeeded, while final consensus was correct for
only 16 pairs. Other relationship classes remain untested by this diagnostic.

Cross-review repaired haplotype matching in IBD2 conversion, probability ties,
malformed topology handling and unsupported independent phase blocks. Genotype
preparation rejects unphased heterozygotes and requires an explicit missing-call
policy. New validation records bind protocol, truth, source, model and output
bytes. Historical runs retain their original source snapshots and weaker
provenance contracts; these changes do not retroactively validate repaired
pipelines. The [parallel follow-up](parallel_followup_2026_11_20.md) records the
experiments and the independent review findings.

**Task 12 — independent reference validation.** A fixed corrected refit was
tested on four new pedigree seeds, using calibration-only correction references.
All 40 CSVs passed pair coverage and format checks, but exact accuracy fell from
37.64% to 37.37%. All four seed-average changes were negative. The candidate is
not promoted; the component experiment does not certify joint family
reconstruction. A later arithmetic review found that the legacy correction
subtracts short segments from only the twenty longest stored lengths, corrupting
the full histogram. A separate versioned repair preserves old model semantics;
the fixed refit and four new confirmation seeds have now completed. Exact
accuracy is 37.14%, versus 37.51% for the historical comparator and 37.47% for
the previous refit. Both prespecified relative gates fail. All 28,000 template
targets remain in the full-template denominator; 27,260 have supplied segments.
The reported accuracy above is conditional on supplied segments. The arithmetic
repair is retained, but the candidate is not promoted. See the [Task 12 report](task12_readiness_2026_11_20.md) and
[histogram repair review](task12_histogram_repair_2026_11_20.md).

## Scientific limits

Task 5's earlier 25.4 kb median localization error was retrospective and is not
fresh release evidence. Tasks 6 and 8 fail important localization/detection
checks. Tasks 1–4 require broader model recovery, artifact regeneration where
provenance was invalidated, and uncertainty validation. Task 7 lacks matched
BGS validation. The complete Tasks 10 and 11 diagnostics expose accuracy
failures; repaired pipelines need new validation. Task 12 lacks complete matched
pipeline validation. Pair-template performance is insufficient to certify
family reconstruction.

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

The initial readiness repair passed 340 tests and was merged in
[PR #4](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/4).
The combined source in [PR #5](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/5)
passes all **949 tests** on macOS/Python 3.12,
Linux ARM/Python 3.12 and Betty/Linux x86-64/Python 3.11. Ruff and GitHub CI also
pass. It integrates the completed submission and pedigree reviews; later
Tasks 10–11 and 13–14 development is separate from that verified snapshot.
A separate Task 7/12 repair snapshot passes **1,022 tests** on the same three
platforms, with Ruff passing. This includes two new regressions for a
floating-point tie in the Task 12 comparison scorer, corrected using integer
correct-count differences. Exact source checks and test logs are recorded in the
[repair verification](../experiments/readiness_2026_11_20/repair_integration/verification.json).
Task 1 production optimizer changes and other-chat development are excluded
from that snapshot. Cross-review
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
