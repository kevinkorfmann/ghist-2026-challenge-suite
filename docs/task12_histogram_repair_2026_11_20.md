# Task 12 histogram correction review

Task 12 is not ready for a complete-pipeline reliability claim. A separate,
versioned correction fixes a reproducible histogram defect. The fixed refit and
fresh component confirmation protocol is documented below; no candidate is
automatically promoted. Historical artifacts and prior confirmation results
remain unchanged.

## Confirmed defect and repair

`correct_records` passes the global count offset to
`trim_shortest_background_shape`, which has only the twenty longest stored
segments. It removes the shortest entries *within that truncated list* and then
subtracts their bin counts from the complete histogram. Those are not generally
the shortest segments in the input.

For a clean 60 cM segment plus one hundred 6 cM background segments, correction
against the clean reference correctly reduces the marginal count to one and
total length to 60 cM. However, the old histogram retains eighty-one 6 cM
segments and deletes the 60 cM segment. The Task 12 model excludes ordered-length
features at its 5 cM threshold, but it still uses these corrupted histogram bins.
Training and evaluation share this transform; its contribution to measured
accuracy cannot be inferred from this arithmetic reproduction alone.

The new `scripts/task12_background_correction.py` subtracts count mass from the
complete histogram, starting at the shortest eligible bin. Fractional quantile
offsets remain fractional; exhaustion saturates at zero. It validates complete,
finite, nonnegative raw bin/count/length features and input consistency. The
existing separate quantile corrections to counts and total lengths are retained.
Their resulting marginal summaries are estimates, not a reconstructed physical
segment list or a claim of consistency at every threshold.

The API deliberately distinguishes `task12-top20-v1` from
`task12-full-histogram-v2`. `correct_for_model` reads artifact metadata under
`task12_background_correction`, containing `version`, `minimum_cm`, and
`genome_cm`. Missing metadata selects only historical v1; malformed or unknown
explicit metadata raises an error. V2 requires exactly the existing 5 cM
histogram/count/length feature list and weight 1.0. Ordered/max-length features
and blending a historical baseline are rejected because their semantics are not
covered by this repair. Existing training and prediction entry points are not
rewired. A v2 model must be explicitly fitted and marked before using the API.

The initial correction-only review passed 53 focused tests. Coverage includes the >20-segment regression,
400 comparisons against complete-length histogram oracles, fractional and
saturated removal, subthreshold preservation, input immutability, invalid domains,
and metadata/feature/blend compatibility. Ruff passed. The subsequent bounded experiment adds Task 12-only runners and provenance/
scoring tests. Task 9, Task 10, Task 11 and shared legacy correction sources are
unchanged by this work.

## What the preserved failure evidence shows

The prior four-seed confirmation report SHA-256 is
`fd9e8ad64ef18d492833e75e5d8ebc235fde7de4a120c886b53a60fdfbf22fd1`.
All forty stored CSV hashes were checked before read-only confusion summaries.
The independent-reference candidate's equal-seed/regime mean class recalls are:

| Relationship | Recall |
|---|---:|
| Siblings | 99.95% |
| Avuncular | 97.61% |
| Avuncular once removed | 50.71% |
| First cousin | 53.33% |
| Avuncular twice removed | 37.89% |
| First cousin once removed | 37.20% |
| First cousin twice removed | 21.20% |
| Second cousin | 21.24% |
| First cousin thrice removed | 10.34% |
| Second cousin once removed | 8.01% |
| Second cousin twice removed | 11.54% |
| Third cousin | 7.14% |
| Third cousin once removed | 19.44% |
| Fourth cousin | 47.62% |

In pooled repeated-regime outputs, 39.19% of avuncular-once-removed pairs are
predicted as first cousins, and 36.04% of first cousins as avuncular once removed.
The corresponding confusion between avuncular twice removed and first cousins
once removed is 32.45% and 32.90%. These labels have the same meiosis count and
expected sharing, so aggregate sharing alone cannot distinguish them; their
distributions and family context can differ. This is not evidence of a label-map
swap. Cross-meiosis errors still account for 46.47 percentage points, so this
ambiguity does not explain the whole failure.

Earlier matched-seed historical-model results were already only 47.52% exact
accuracy without added contamination, versus 38.49% at 0.6× and 36.36% at 1.4×
contamination. That is a different fixed model/seed comparison, not a new clean
evaluation of the refit. The poor clean baseline limits what a background-only
repair can be expected to establish.

The four confirmation definitions generated 28,000 target pairs, of which
27,212 appear in exact-segment files. Missing zero-IBD pairs include 527 of 2,000
fourth cousins (26.35%), 194 third cousins once removed, 33 third cousins and
31 second cousins twice removed. The other three omitted pairs are one first
cousin thrice removed and two second cousins once removed. CSV coverage is
correct for the supplied segments; it is not complete template or VCF coverage.

The balanced isolated-pair template mixture, global 25th-percentile correction,
independent added background segments, and inherited unlabeled empirical profile
do not establish transfer to a full family with different relationship frequencies
and correlated founder ancestry. Exact segments also bypass detector, genotype,
phasing and topology errors. No new empirical or challenge inputs were read.

## Fixed bounded experiment

The original [proposed protocol](../experiments/readiness_2026_11_20/task12_histogram_v2/proposed_protocol.json)
remains preserved with SHA-256
`e61cccf001e09b24a6d8c045a1cab741057f31276b8d8568ec5888d076412759`.
Its executable successor is [frozen_protocol.json](../experiments/readiness_2026_11_20/task12_histogram_v2/frozen_protocol.json),
SHA-256 `ce2a51405f323807bd9c8e1eebe20ce7ea706b20e7b3176ef2e6b774ca3ad5ed`.
The source and input hashes include 23 selected code, test, configuration and
audit files, the original six calibration simulations and their seed logs,
two historical comparator models, the baseline model, and simulator/map inputs.
The separate runtime manifest hashes the Python executable and 3,630 installed
files from six scientific/runtime packages. It does not claim to hash every
standard-library or operating-system library.

The single fixed fit uses the original 240-iteration HGB, blend 1.0, and six
calibration replicates. Each replicate's correction reference consists only of
the other five replicates. Architecture, features, augmentation settings, class
weighting and random state remain fixed. The runner requires exactly 203,950
training rows and the original 24 augmentation settings. Calibration checks
validate representation and provenance; they do not select a model by training
accuracy. Candidate metadata explicitly marks v2 correction and new training
provenance; arbitrary historical metadata and validation scores are not inherited.

Only after a candidate artifact, training record and lock exist may the runner
generate four new Ped-sim seeds, 2026115011–2026115014. Calibration seed logs
101226000–101226005 are hash-checked and parsed at startup; the executable
registry enforces distinct calibration, confirmation and augmentation streams.
The 24 new simulation/augmentation seeds had no matches in the bounded audit
of 699 local and 2,915 Betty source/protocol/log files. Exact scope and exclusions
are preserved in the audit files. This is not a claim about undisclosed future
plans in other chats. The prior four confirmation seeds are excluded.

All five original regimes are fixed. Each new augmented input is evaluated with
both the historical model and the previous independent-reference refit under
their original v1 correction, as well as the candidate under v2. Reports retain
the generated 7,000 target pairs per seed separately from segment-present pairs.
Missing zero-IBD predictions count as failures in the complete-template metric.
Regimes repeat measurements on the same four simulated pedigrees; they are not
twenty independent confirmation replicates. CSVs pass the existing Task 12
validator with exact supplied-pair coverage.

The final focused suite passed 74 tests in 14.35 seconds; Ruff passed. Added
regressions cover a protocol change during the real training entry point,
calibration seed-log tampering/reuse, candidate lock tampering, exact probability
and pair coverage, and an emitted non-first tied maximum. The tie change
generalizes the scorer contract; it is not evidence of a bug in current model
outputs. Both runners bind parsed protocol bytes to a captured digest and reject
changes before/after fitting, generation and final evidence publication.

Three earlier jobs are preserved as cancelled attempts: 10851898 stopped after
two seconds, training 10851952 after forty seconds, and its dependent confirmation
10851959 before starting. Final inspection of 10851952 found only an empty
candidate directory, preserved in the superseded snapshot. No fitted model,
training report, candidate lock or confirmation data exists from those attempts.
Because fit entry was not logged, the empty directory does not prove that fitting
was never invoked. An earlier immediate cancellation check saw no directory;
[cancellation_stage_correction.json](../experiments/readiness_2026_11_20/task12_histogram_v2/cancellation_stage_correction.json)
records this correction and preserves the superseded metadata. No result-based
model or seed changes occurred.

Training job 10852084 completed with exit 0 in 1m53s (peak RSS 3.21 GB).
Confirmation job 10852118 completed with exit 0 in 2m23s (peak RSS 4.76 GB).
Both respected separate fifteen-minute, eight-CPU Slurm limits. The model and
training/lock evidence were fixed before any fresh confirmation simulation.

## Fresh confirmation result: candidate rejected

The histogram repair did not improve predictive accuracy. Equal-seed/regime
averages from the preserved primary report are:

| Metric | Historical model | Prior reference refit | Histogram-v2 candidate |
|---|---:|---:|---:|
| Exact label accuracy, segment-present pairs | 37.5133% | 37.4699% | 37.1416% |
| Balanced accuracy, segment-present pairs | 37.5973% | 37.4853% | 37.2055% |
| Meiosis accuracy, segment-present pairs | 53.3903% | 53.3698% | 53.1507% |
| Full-template accuracy, omitted pairs wrong | 36.5221% | 36.4800% | 36.1600% |
| Negative log likelihood | 1.494355 | 1.500824 | 1.508542 |
| Multiclass Brier sum | 0.669183 | 0.671421 | 0.674326 |
| Ten-bin ECE | 0.026207 | 0.032783 | 0.039611 |

Both prespecified relative gates failed. Exact accuracy declined by 0.3717
percentage points versus historical and 0.3283 points versus the prior refit.
Against historical, three seed means declined and one tied exactly; against
the prior refit, three declined and one improved. The worst individual
seed/regime exact decline was 1.4661 and 1.1582 percentage points respectively.
Regimes are paired repetitions within four independent simulations, so these
are descriptive results, not a twenty-replicate confidence claim.

The four seed denominators are 6,833, 6,797, 6,821 and 6,809 segment-present pairs,
27,260 of 28,000 generated targets. The 740 omissions comprise 511 fourth
cousins, 180 third cousins once removed, 27 third cousins, 19 second cousins
twice removed, two first cousins thrice removed and one second cousin once
removed. All 60 prediction CSVs passed production schema/pair validation.
A separate local audit verified their hashes, rederived every confusion matrix
from the simulated pair labels, and checked identical pair sets across all
models/regimes for each seed. The candidate remains unpromoted.

The new candidate's recall is still 99.89% for siblings and 98.16% for avuncular
pairs, versus approximately 7.55–22.09% for several distant-relationship labels.
The arithmetic defect is real, but this experiment does not establish it as the
cause of the low overall accuracy. Histogram correctness alone does not resolve
relationship ambiguity, background-model mismatch, or absent family context.
The four fresh seeds are now observed diagnostic evidence and cannot certify
another revised model. No further tuning or model experiment occurred.

## Preserved primary evidence and separate scorer correction

The primary confirmation SHA-256 is
`78b1a925ddb2dc710a6c23089443b08cbffeb303299f9d1fc90f61b1e900fa64`.
The candidate model SHA-256 is
`de2d3aa1a96fc778aced63ea10754894c4990f38473496fcccd91540fbf3ab48`.
It is preserved on Betty under
`reviews/readiness-task12-histogram-v2-20261120/artifacts/task12_histogram_v2_20261120/`
and backed up locally under
`/Users/kevinkorfmann/.cache/ghist-2026/task12-histogram-v2-20261120/`.

The primary scorer's floating subtraction/averaging reported a gain of
`1.1102230246251566e-17` for historical comparison seed 2026115012 despite
identical integer correct totals. The five candidate totals were
2632, 2580, 2478, 2433, 2459; historical totals were
2617, 2546, 2496, 2447, 2476, each using 6,797 pairs per regime.
The prospective scorer now compares integer diagonal-count differences before
dividing, and rejects unequal paired/regime denominators. It correctly treats
that seed as a tie. Two new regression tests bring the focused suite to 76;
Ruff passed, and root's local full integration passed 1,022 tests.

The original report/source was not rewritten. The 23 frozen source/configuration/test/audit
files plus protocol/runtime/training/lock are preserved in
[primary_run_source.tar.gz](../experiments/readiness_2026_11_20/task12_histogram_v2/primary_run_source.tar.gz),
SHA-256 `179e652a4c099028d13793a5f07c0b16660e7e07493775e6506dae5a8034b0c0`.
The source hashes inside `frozen_protocol.json` intentionally describe that
primary runtime, not the prospectively corrected scorer in the repository.
The Betty primary snapshot remains unchanged. All 60 CSVs and 20 per-regime
checkpoints are preserved in
[prediction_csvs_and_checkpoints.tar.gz](../experiments/readiness_2026_11_20/task12_histogram_v2/prediction_csvs_and_checkpoints.tar.gz),
SHA-256 `26d327ec29afaff1f24d290f14bcbae7be9e6d19ff96371f4630dd19cb1683c6`.

[comparison_recomputation.json](../experiments/readiness_2026_11_20/task12_histogram_v2/comparison_recomputation.json),
SHA-256 `7e4d4d639885db0f8e028097d830e67aa464440a8d9830f9dc69e73179c41f57`,
records both old and corrected comparisons, exact source/report identities and
all verified CSV hashes. Every original overall gate decision and aggregate
metric is unchanged. Its companion `recompute_comparisons.py` reproduces the
separate arithmetic audit from the frozen primary report/protocol and extracted
CSV archive. No model was refitted, no predictions regenerated, and no candidate
selected after viewing these outcomes.

This remains exact-segment component validation. The full founder-VCF,
IBD-detection and pedigree/topology route remains unvalidated. No current
challenge inputs, public answers or competition submissions were used.
