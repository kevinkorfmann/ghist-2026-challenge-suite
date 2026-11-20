# Task 11 consensus correctness audit, November 20

**Task 11 remains unready.** Deterministic bookkeeping and input-validation defects
are repaired, and the historical freeze now requires explicit independent
complete-task evidence. No new prediction policy has been selected for production,
no final challenge inputs were used, and no submission was made.

The original failure remains **16/96 correct after consensus versus 96/96 before
consensus at every HBD threshold**, across 16 sibling-only synthetic families in
two cohorts. See the unchanged
[whole-pipeline report](task11_complete_pipeline_2026_11_20.md). The two cohorts
were already inspected when this audit began and are now development data for
any policy designed from these findings.

## Deterministic repairs

[`consensus_task11_dog_pedigrees.py`](../scripts/consensus_task11_dog_pedigrees.py)
now rejects duplicate unordered pairs, self-pairs, inconsistent HBD pair sets,
duplicate family samples/IDs, incomplete within-family coverage, missing topology
pairs, missing/unknown probability labels, invalid/nonfinite probabilities and
unnormalized distributions. Finite fractional topology distances can no longer
be silently truncated into an accepted relationship. Nonsequential family IDs no
longer disable sibling-component propagation.

Every reported `confidence` is now the **mean original forest probability for the
final label**, including after relation-vector refinement and propagation. It is
explicitly uncalibrated; it is not a probability that the complete inference is
correct. Constrained probabilities are no longer reported as though they were
unconditional confidence. A constraint with exactly zero probability mass raises
an error instead of selecting an arbitrary zero-support label.

`hbd_sensitivity_stability` and per-row `sensitivity_labels` now describe the raw
forest argmax before topology replacement. The old conditioned statistic is
reported separately as `topology_conditioned_hbd_sensitivity_stability`. The
historical refinement rule still uses conditioned labels; this reporting repair
does not silently replace its inference policy. An empty fallback subset yields
JSON `null` with a zero denominator, rather than nonstandard `NaN`.

The old topology overrides, meiosis constraint and propagation remain in the
legacy method. Their accuracy has not been repaired by these bookkeeping changes.
The replay below used `task11-legacy-consensus-audited-v1`, development only.
The subsequent independent cross-review advances the current source to
`task11-legacy-consensus-audited-v2`; it does not rewrite those replay results.

## Replay and the new guard failure

[Replay evidence](../experiments/consensus_audit_2026_11_20/task11/audit_summary.json)
and [confidence corrections](../experiments/consensus_audit_2026_11_20/task11/confidence_corrections.csv)
preserve the original results. All **42 retained outputs with original recorded
hashes** verify, and all **49 historical local files** remain byte-for-byte
unchanged.

| Development cohort | Audited legacy behavior | Guarded provisional pairs | Abstentions |
|---|---|---:|---:|
| 2026112711 | Runs; all 48 labels unchanged, 7/48 correct | 3/48 | 45/48 |
| 2026112811 | Refuses a zero-support topology constraint | 4/48 | 44/48 |

The refused constraint is HBD10 for `synthetic11_22` / `synthetic11_23` in the
second cohort. Bonsai supplies `[0, 3, 1]`, but the forest assigns probability
1.0 to siblings and zero to **every** three-meiosis label. The historical
sensitivity code nevertheless selected one. The audited legacy script now
produces no candidate for this cohort; the audit retains the original labels
and computes their probability support separately, clearly marked as metadata
analysis of historical output.

Across both historical outputs, three propagated `half_sibling` calls had stale
confidence 1.0 despite final-label forest support 0. Three other calls reported
conditional confidence 0.809–1.0 while their unconditional support was only
0.000923–0.00722. The metadata audit records all six corrections. Raw HBD stability
is still 1.0 on these particular cohorts because the forest was stable; a focused
test with differing raw calls proves that a common topology override no longer
hides that difference in the reported raw statistic.

## Separate development policy and tradeoff

[`consensus_task11_guarded_candidate.py`](../scripts/consensus_task11_guarded_candidate.py)
requires explicit selection of
`task11-guarded-pair-proposal-development-v1`. It preserves the unique argmax of
the averaged original forest probabilities as a proposal. A relationship is
provisional only when at least two HBD runs have the same unique argmax and a
mapped topology label agrees. Ties, HBD disagreement, unassessed sensitivity,
unmapped topology or topology conflict cause abstention. Topology cannot replace,
constrain or propagate a proposal. The tool writes development JSON only and
always reports `submission_ready=false`; it is not selected by the default or
freeze path. Provenance binds its own source, the imported consensus helpers,
the label definitions and all four input files.

The [fixed development protocol](../experiments/consensus_audit_2026_11_20/task11/policy_protocol.json)
was written before this replay. No probabilities were fitted and no numerical
thresholds were tuned. All 96 forest proposals match the sibling-only truth
retrospectively, but only **7/96** pairs remain provisional; **89/96 abstain**.
Counting abstentions as failures, the guarded policy therefore yields **7/96**,
not 100% whole-task accuracy. Agreement is not evidence of calibration.

Using topology solely as advisory information and exporting every forest
proposal would improve coverage on these examples, but would select an
uncertified method from the same examples used to design it. Abstention instead
exposes the unresolved model conflict. Its severe coverage loss makes this a
development diagnostic, not a usable submission method.

Before promotion, a prospectively frozen confirmation protocol must evaluate
new independent whole pedigrees and founders, all eight labels, varied
topologies/inbreeding, matched dog founder LD/frequencies and permitted
missingness. It must score all expected pairs, abstentions and missing pairs as
failures, spurious cross-family pairs, per-label confusion and complete-pipeline
balanced accuracy, with family/cohort-level uncertainty. Acceptance criteria must
be recorded before those new outcomes are examined. More sibling-only replicates
cannot establish whole-task performance.

## Freeze guard and provenance

[`freeze_task11_multifamily_inbreeding.py`](../scripts/freeze_task11_multifamily_inbreeding.py)
now stops before historical artifacts, genotype inputs or writes unless
`--validation-evidence` is supplied. The `task11-complete-validation-v1` evidence
must identify the selected method, all eight evaluated labels and new confirmation
cohorts; attest to independent full-pipeline validation, development separation,
varied topology, matched founders/missingness, unconditional scoring and passed
predeclared criteria; and bind protocol/report, candidate, summaries, all three
models/prediction sets, relevant source, map, hap-IBD JAR and inference configuration
by SHA-256. The known development seeds are explicitly ineligible for confirmation.
The source list pins feature extraction in `fit_task11_dog_relationships.py`, the
shared conversion helpers and relatedness modules; there is no separate
`features.py` module in this pipeline.

The guard checks **provenance plus explicit scientific attestations**. It cannot
machine-certify the truth of those attestations. No qualifying evidence was
created during this audit. Once evidence exists, the existing robustness checks
still run, now using raw HBD sensitivity, an explicit empty-subset rule, the
audited method version and nonzero final-label probability support. Guarded
development output cannot pass by being silently substituted for the selected
legacy method.

## Adapter finding and limits of attribution

The hash-verified
[original adapter snapshot](../experiments/consensus_audit_2026_11_20/task11/original_source/scripts/reconstruct_bonsai_pedigree.py)
matches the smoke protocol's SHA-256
`5d317416d437a92746af4c618001af2d521b2854d55accba1ac62142bcc4ea4f`.
Its `read_haplotype_segments` discarded haplotype IDs at original lines 41–47;
`coverage_intervals` marked any multiplicity of at least two as full sharing at
line 62. Two matches sharing one haplotype do not establish two disjoint diploid
matches, so this can manufacture IBD2 in inbred data. The Task 11 adapter also
scales dog map lengths by `3544.57 / dog_genome_cm` before calling Bonsai; that
scaling alone does not establish that the retained likelihood model is appropriate.
These are plausible contributors alongside the unconditional overrides, not an
established single cause of the observed topology errors.

A separate agent is repairing conversion/feature extraction. This audit used
the original stored forest probabilities and topology throughout. Its replays
therefore cannot certify the changed current genotype-to-feature pipeline, and
old fitted-model performance cannot be assumed to survive that feature change.

Validation: **54 focused tests passed** across consensus behavior, freeze guards
and the existing whole-pipeline smoke harness; Ruff passed for the changed code
and replay script. Tests use synthetic fixtures, not purported release evidence.
Reproduce the audit with
`PYTHONPATH=src:. .venv/bin/python experiments/consensus_audit_2026_11_20/task11/replay_audit.py`.

## Independent source cross-review after the recorded runs

The follow-up review reproduced an ordering defect: identical probability
vectors with a 0.5/0.5 sibling/half-sibling tie reported HBD stability 1.0 with
matching JSON key order and 0.0 with reversed keys in one run. The current v2
source resolves exact ties using the declared `MEIOSES` label order. This is
deterministic bookkeeping, not stronger biological evidence; the separate
guarded development policy continues to abstain on ties.

The review also reproduced unconstrained fallback from malformed tuples such
as `[NaN, 1, 2]`. Current source rejects these and permits only the explicit
`(+Infinity, +Infinity, None)` unrelated sentinel used by
[Bonsai's pedigree object](https://github.com/23andMe/bonsaitree/blob/master/bonsaitree/pedigree_object.py).
The standard mapped finite relationships are unchanged.

The combined follow-up suite passed **170 tests**, including consensus, freeze,
genotype phase contracts, segment conversion, missingness and both synthetic
complete-pipeline harnesses. An independent discrete-coordinate oracle matched
the corrected segment converter in **5,000 randomized interval sets**. Ruff and
diff checks passed. All **171 preexisting experiment files** in the complete-
pipeline and consensus-audit evidence directories retained their hashes.
[Separate cross-review verification](../experiments/complete_pipeline_2026_11_20/cross_review/verification.json)
records the before/after source identities and the preserved-artifact manifest.
No full-pipeline simulation or inference was rerun with v2, and no accuracy
improvement or release readiness is inferred from these software checks.
