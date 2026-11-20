# Readiness follow-up, 20 November 2026

Task 9 now has an explicitly selected, versioned release candidate that passes
its scientific and production checks in the tested exact-IBD input domain.
This is the strongest working challenge pipeline. Tasks 1–8 and 10–12 still
have scientific validation or accuracy gaps. The two later low-coverage
challenges are being developed in a separate, unmerged continuation. No current
final-data inference or competition submission has been performed by this review;
separate submission work is outside its frozen validation scope.

The user authorized direct merging and parallel scientific review, while asking
that final-data runs wait. Three reviewers cover demography, sweeps and
relatedness; the coordinating reviewer covers Task 1 and integration. This
report supplements the earlier [scientific review](scientific_validation_2026_11_20.md),
whose archived results and historical artifacts remain preserved.

## What is working

| Challenge | Current evidence | Status |
|---|---|---|
| 1 — bottleneck | Export checks pass; versioned optimizer repair passes 1,089 tests on three platforms; original 24-case panel has major recovery errors | Not ready; recovery and uncertainty gaps |
| 2 — secondary contact | All original cases finish; malformed interval exports are rejected; hidden-history and contact-time errors remain | Runs, not scientifically ready |
| 3 — growth/BGS | New 10 Mb neutral simulations recover strong growth; weak growth remains poorly identified | Runs, broader validation needed |
| 4 — split/migration/BGS | Numerical stall fixed; Q1 fit finishes, but split-time error is +47% | Runs, accuracy/calibration gaps |
| 5 — hard sweep | Three completed cases: 2/3 interval coverage, 0/3 agreement passes; default-policy retry exhausted 22 GiB and both fixed-interval trials timed out | Not ready; simulation feasibility and recovery gaps |
| 6 — soft sweep | Eight new development configurations fail recovery criteria | Not ready |
| 7 — sweep/BGS | Canonical export, real RAiSD, prediction and validators pass on one retained development tree; matching biological validation missing | Not ready |
| 8 — multiple sweeps/BGS | Existing detection/localization failures remain | Not ready |
| 9 — exact IBD relationships | Fixed candidate passes scientific/production gates; separate submission has final weighted-Jaccard scores 0.834/0.842 | Validated candidate in tested domain; submitted |
| 10 — human pedigrees | Complete pipeline produces all 153 pairs in one family; only 58 exact labels correct | Runs, accuracy inadequate |
| 11 — dog pedigrees | Two cohorts: classifier 96/96 sibling pairs correct, final consensus 16/96 | Consensus fails; broader relationships unvalidated |
| 12 — background-related families | Full-histogram repair passes software checks but fails fresh accuracy criteria | Not ready |
| 13 — low-coverage bottleneck | Separate continuation reports a full-width BAM run, with optimizer, target-stability and fraction-recovery failures | Development, not ready |
| 14 — low-coverage split | Development passed, but first fresh confirmation cohort failed seed stability; separate stabilization experiment registered | Current policy fails confirmation |

The Tasks 13–14 rows summarize separate continuation evidence. That unmerged
implementation is outside the frozen integration snapshots below. Task 14’s
first confirmation cohort changed one target by 1.2769% across estimator seeds,
exceeding the unchanged 1% stability limit. This prevents passage of the original
all-six policy; its separate stabilization experiment does not erase that failure.
The [external-work observation](../experiments/readiness_2026_11_20/followup_integration/continuation_observation.json)
binds the reports used for this status, without claiming a new model evaluation.

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
input and prediction fields explicitly pending in that historical offline handoff.
Separately, submission 9785353 now has final scores of **0.834 for relation** and
**0.842 for meioses count** on the [official dashboard](https://www.synapse.org/Synapse:syn74803798/wiki/640981).
These are truth-class-size-weighted Jaccard scores, not ordinary accuracy.
The separate submission record and its dashboard evidence hashes were checked;
these results did not change the frozen model or its qualification.

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

**Task 2 — interval export integrity.** Independent actual-entrypoint
reproductions showed that empty, reversed and negative-population intervals
could pass the secondary-contact freeze guard. The repair requires the exact
five producer fields, two finite ordered endpoints, positive population/split
values, and nonnegative contact-time/migration values. Valid zero boundaries
remain accepted. Four valid independent checks, including producer-generated
intervals, preserve the original export bytes; failure checks preserve the
existing output and leave no temporary artifacts. The
[repair record](../experiments/readiness_2026_11_20/secondary_contact_interval_integrity/review.md)
and [independent review](../experiments/readiness_2026_11_20/secondary_contact_interval_integrity/root_independent_review.json)
describe the bounded scope. These shape/domain guards do not establish the
planned bootstrap denominator, empirical coverage, or biological recovery.

**Tasks 3–4 — prediction and bootstrap integrity.** The freeze wrappers now
recompute projection medians from every required converged fit, reject malformed
or nonfinite intervals, and validate the available theta/scaled-to-physical
arithmetic for each observed projection and the refined Task 4 fit. Componentwise
medians remain aggregates and are not treated as individual model fits.
Bootstrap producers declare the complete planned denominator before fitting;
summaries require every unique planned replica, a converged selected optimum,
consistent inputs/source/settings, and retained fit arithmetic. Freeze recomputes
the interval statistics from those records. Historical summaries lacking this
contract must be regenerated; failed replicas cannot be omitted to authorize an
interval. Strict atomic output preserves existing artifacts on failure.

Independent reviews reproduced two additional consistency gaps in intermediate
patches: corrupted physical timing in retained bootstrap fits, then corrupted
observed projection timing with a matching cached median. Both are fixed in the
final source; the earlier passing test runs and failure reproductions remain
preserved. These checks establish artifact consistency, not interval coverage,
full selected-endpoint provenance, or recovery accuracy. See the
[integrity report](demography_freeze_integrity_2026_11_20.md) and its
[final independent review](../experiments/readiness_2026_11_20/demography_freeze/root_independent_reviews/observed_v3/review.json).

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

All 24 cases and 36 method problems completed under the original frozen coarse
optimizer. Only 30 selected fits report convergence: 8/12 matched moments,
10/12 matched dadi and 12/12 ancestral-history moments. Matched median absolute
relative errors are 7.35% and 8.82% for time, 1.24% and 1.33% for the decline
fraction, and 0.69% and 0.84% for ancestral size (moments and dadi respectively).
These medians conceal the weak-decline case 16: time errors are +725% and +673%,
and ancestral-size errors are +2,209% and +1,424%, despite reported convergence.
Moments also flags poor restart agreement. Under the perturbed ancestral
history, median time error is 21.16%; case 23 reaches +2,884% time and +4,675%
ancestral-size error. These are diagnostic results on six parameter points,
not a full-prior accuracy estimate. All unfavorable outcomes remain in the
[completed original panel](../experiments/readiness_2026_11_20/bottleneck/full_panel_verification.json).

The original VCF identity records remain rejected as legacy records; the
separate corrected v2 export proof above supplies the verified export evidence.
The original scientific results are not relabeled as outputs from the repaired
optimizer. The [production optimizer review](bottleneck_optimizer_production_2026_11_20.md)
records the coarse/fine/refinement development, source replay, residual agreement
failure, and separate stress retry. Numerical repair alone does not establish
identifiability, reliable biological recovery, or interval coverage.

After the final source and all replay outcomes were fixed, a separate
[descriptive truth comparison](../experiments/readiness_2026_11_20/bottleneck/production_stress_biological_context.json)
confirmed that the repaired methods still misestimate the known stress case:
about 218,057 and 216,381 generations, versus the simulated 28,000, despite
passing the numerical gates. This single case was selected after its earlier
poor recovery was known; it is development evidence, not fresh confirmation.
No model, setting, or acceptance threshold was changed after this comparison.

A subsequently frozen [four-point diagnostic](bottleneck_stress_fixed_points_2026_11_20.md)
evaluated only the generating and selected shapes for each engine on that same
existing spectrum. Selected likelihoods, theta and physical transforms reproduce
exactly. Selected shapes score 9.802 and 7.556 likelihood units above the true
shapes, while their current-size products are within 0.56–0.59% of truth. This
excludes a physical-transform discrepancy in the checked values; it does not
explain the likelihood preference, establish global optimality or identifiability,
or provide independent recovery evidence. No optimizer or parameter search ran.

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
A separate paired synthetic software fixture checked more frequent ancestry
simplification as a possible resource remedy. Both arms retained the same
population/allele trajectory, sampled identities, selected genotypes and
canonical marginal ancestry, including a forced checkpoint restore. Its original
initialization-syntax failure and corrected run are both retained. This tiny
fixture provides no credible full-case memory-savings estimate or biological
recovery evidence; production defaults remain unchanged. See the
[fixture review](../experiments/readiness_2026_11_20/sweeps/task5_simplification_fixture_review.json).
A separately frozen full-case attempt then used fixed interval 20 under the same
4-CPU, 22-GiB, 45-minute bound, with the original case 5 seed, chromosome,
demography, maps, scaling and acceptance conditions unchanged. Slurm job
10855907 ended in **TIMEOUT** after 45 minutes 9 seconds. Its 270 resource samples
record a 2.023 GiB high-water mark, no OOM or allocation-failure counters, and a
maximum sampling gap of 10.002 seconds. The final sample precedes scheduler end
by about 7.14 seconds; termination prevented the adapter's finalization, leaving
a stale `RUNNING` record that does not indicate success. The last logged tick is
27,200 with 43,410 diploids, still in the first growth phase; no checkpoint
save/restore, final simulation, recapitation, neutral overlay, VCF export or
scanner stage completed. Later population peaks of about 71,600 remain
unobserved. The [terminal collection](../experiments/readiness_2026_11_20/sweeps/task5_interval20_feasibility/terminal/terminal_collection.json)
binds the raw traces and checks all 22 frozen source files and 148 external
locks. This censored run provides neither a full-case memory upper bound nor a
matched-phase memory-saving estimate. The original OOM and this timeout both
remain failures. A [prospective independent review](../experiments/readiness_2026_11_20/sweeps/root_independent_reviews/task5_interval100_prospective/review.json)
sets explicit conditions for a separate interval-100 software fixture and
bounded resource pilot; it is not evidence that either succeeds.

**Task 5 — interval-100 resource trial.** The separately registered tiny
20-versus-100 simplification fixture passed its trajectory, sample-identity,
selected-genotype and canonical-ancestry checks, including a forced restore.
Independent checking confirmed 2,208 sample-pair ancestry comparisons. It is
software evidence from one reused seed, with no full-case performance claim.
The subsequent same-seed full-case interval-100 trial, job 10857763, also ended
in **TIMEOUT**, after 45 minutes 10 seconds under the unchanged 4-CPU/22-GiB
bound. Its 271 resource samples record a 7.791 GiB high-water mark with zero
recorded allocation failures or OOM kills. The largest sample gap is 10.002
seconds; the final 1.540 seconds before scheduler end are unsampled.
It reached the first contraction at tick 33,600 with 52,855 diploids, after a
maximum logged population of 71,194. It did not finish SLiM, save/restore a
checkpoint, produce a final tree, or reach recapitation, mutation overlay,
VCF export or scanning. The retained adapter `RUNNING` status is stale.
All 24 frozen source files, 148 external locks, runtime and midpoint trace
prefixes passed checks. The [terminal collection](../experiments/readiness_2026_11_20/sweeps/task5_interval100_feasibility/terminal/terminal_collection.json)
records the partial phase coverage. This run reaches a later phase than the
interval-20 attempt but establishes neither a full-case memory bound nor a
matched-phase performance improvement. Both timeouts and the earlier OOMs
remain failures; no new biological replicate or accuracy evaluation is added.
Production defaults and scientific acceptance thresholds are unchanged.

The tiny completed subset cannot establish full-panel recovery. See the
[sweep report](readiness_sweeps_2026_11_20.md) for failures and resource policy.

**Task 7 — actual retained-tree pipeline replay.** One previously generated
development tree now passes canonical VCF export, real RAiSD scanning, production
prediction, BED validation and candidate recomputation. The 25-diploid replay
retains 777 exported sites and 728 callable RAiSD windows. The original attempt
failed before export because tskit refuses to trim a tree with retained reference
metadata. A separately reviewed harness retry clears only that unused metadata
on an in-memory copy, verifies unchanged biological tables and genotypes, and
checks the coordinate shift. The raw tree and original failure remain preserved;
all eleven production source files are unchanged. The retry completed in nine
scheduler seconds. The 1 Mb segment is too short for a reference peak more than
2 Mb away, so the expected consensus rejection and absent signal ratio remain
intact. This is one synthetic-segment interoperability check using an existing
stress tree, not a new independent simulation, matched-BGS recovery estimate,
whole-chromosome accessibility check, freeze or release qualification. The tree
digest is a retrospective attestation, not an original pre-simulation lock.
See the [terminal report](../experiments/readiness_2026_11_20/task7_development_replay/retry_v2/terminal_report.md).

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

Claude access recovered at 21:00 UTC. Two duplicate, read-only source reviews
completed against the frozen Task 1 optimizer/bootstrap and Task 3–4
export/bootstrap code. Their source hashes were unchanged before and after
review; no challenge data, fitted outputs or other reviews were supplied.
The [review records](../experiments/readiness_2026_11_20/claude/final_duplicate_reviews/task1_optimizer_bootstrap.metadata.json)
and [Task 3–4 record](../experiments/readiness_2026_11_20/claude/final_duplicate_reviews/task3_4_export_bootstrap.metadata.json)
preserve the exact scope and outcome. The
[triage](../experiments/readiness_2026_11_20/claude/final_duplicate_reviews/deadline_triage.json)
separates demonstrated candidate-path checks from unverified findings and
deferred interval-validation work. Earlier quota-interrupted attempts still
do not count as completed reviews. These source audits do not establish
scientific accuracy or justify promotion.

The initial readiness repair passed 340 tests and was merged in
[PR #4](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/4).
The combined source in [PR #5](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/5)
passes all **949 tests** on macOS/Python 3.12,
Linux ARM/Python 3.12 and Betty/Linux x86-64/Python 3.11. Ruff and GitHub CI also
pass. It integrates the completed submission and pedigree reviews; later
Tasks 10–11 and 13–14 development is separate from that verified snapshot.
The Task 7/12 repair merged in [PR #6](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/6)
and passes **1,022 tests** on the same three
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

The subsequent Task 1 optimizer and consumer snapshot merged in
[PR #7](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/7) after
GitHub Actions passed, and passes **1,089 tests**
on the same three platforms, with Ruff passing. The final two regressions reject
duplicate original starting vectors and inconsistent paired starts, so repeated
copies of one start cannot create restart agreement. Source hashes for all 465
code/configuration files and 827 ancillary fixture files match on each platform.
The [final integration record](../experiments/readiness_2026_11_20/optimizer_guard_integration/verification.json)
preserves an initial test-snapshot packaging failure: required pre-existing data
fixtures were missing. Restoring those unchanged fixtures resolved the failures;
production source was identical in the failed and successful snapshots. The
prior 1,087-test integration record remains preserved as a separate version.

The Task 3–4 export/bootstrap integrity snapshot merged in
[PR #9](https://github.com/kevinkorfmann/ghist-2026-challenge-suite/pull/9) after
GitHub CI passed and passes **1,260 tests**
with no skips, failures or errors on macOS/Python 3.12, Linux ARM/Python 3.12,
and Betty/Linux x86-64/Python 3.11. Ruff passes on all three. All 468 source/config
files and 827 unchanged ancillary fixtures match before and after each run.
The [integration record](../experiments/readiness_2026_11_20/demography_freeze/integration_v3/verification.json)
links the preserved intermediate snapshots and independent reproductions.
The first Betty harness failed before tests because Ruff was absent; installing
it into an isolated lint directory fixed the harness without changing production
dependencies. Separate Tasks 10–11/13–14 development and the later Task 2 interval
domain repair are outside this snapshot. Scientific readiness is unchanged.

The subsequent Task 2 interval-guard snapshot passes **1,314 tests** with no
skips, errors or failures on the same three platforms; Ruff passes on all three.
All 469 source/configuration files and 827 unchanged ancillary fixtures match
before and after each run. The [final integration record](../experiments/readiness_2026_11_20/secondary_contact_interval_integrity/final_integration/verification.json)
preserves an initial snapshot whose cleanup assertion glob had an extra dot;
the final snapshot corrects only that test line. Independent failure checks used
the correct temporary-file prefix and passed. Neither this integration nor the
Task 7 replay establishes scientific readiness for the remaining challenges.
Separate Tasks 10–11/13–14 development remains outside these frozen snapshots.

Compact protocols, checkpoints, job scripts and verification records are under
[`experiments/readiness_2026_11_20`](../experiments/readiness_2026_11_20/).
Scientific accuracy claims apply only to the stated experiments. No finite test
suite or simulation panel can establish that software is bug-free or that every
inference will be correct.
