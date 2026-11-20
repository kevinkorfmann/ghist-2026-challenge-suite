# Continued readiness work — November 20, 2026

This chat owns new development for Tasks 10–11 and 13–14, plus independent
cohort scoring. The active **Review 2026 ghist challenges** chat continues to
own Tasks 1–8, Task 9 packaging and Task 12. Its shared files, models and
previous experimental results are outside this chat's editing scope.
At 19:13 UTC the user authorized beginning submissions for the strongest
working candidates while other development continues. The coordinator now
owns Task 9 final inference and submission in a separate run directory, using
the completed v2 package unchanged. Other sessions' shared integration remains
outside this chat's editing scope.

Status checked at 20:01 UTC. **Task 9 has been submitted and scored.** Its
current-final-input run completed as job 10856122 in 23 seconds and passed
both production validation and independent identity/coverage review.
Synapse entity `syn75202631` is current version 3, file handle
177273793, MD5 `1edb0a55ced6e19e3a577047185abfa8`, input SHA-256
`ff7031bc00a96f7d38131e010f1168f95fb4a8efad5849c904b2d7e77e4de394`.
The downloaded bytes match the authenticated file sidebar and public entity
metadata. The prediction covers all 997 observed pairs once and was uploaded as
[syn78778815, version 1](https://www.synapse.org/Synapse:syn78778815). It was
submitted as an individual to **GHIST 2026 Relatedness 1 Final**; Synapse displayed
its successful receipt. Submission **9785353** now has official weighted Jaccard
scores **0.834 for relation** and **0.842 for meiosis count**. These are not
ordinary accuracy percentages. The dashboard evidence is retained as
`score_record.json`, `score_dashboard.png` and `score_dashboard.ax.txt`.
Telegram delivery to Kevin’s Saved Messages was visually verified at20:40
Europe/Berlin after the user unlocked the Mac. The output SHA-256 is
`d2810e81b74591c3d10c5305643c3ed4494ab416e44e4f98c37ff6df92b16c6d`.
Run evidence, the independent audit, completed methods record and receipt
screenshot are in `experiments/final_submission_2026_11_20/task9/`.
The other tasks retain the following measured gaps. Rows owned by the other
chat are read-only summaries of its readiness reports, not new experiments.

| Task | Current readiness and remaining gap | Owner |
|---|---|---|
| 1 | Not ready: 240 optimizer endpoints replayed; full-length panel completed with six nonconverged fits and large recovery errors | Other chat |
| 2 | Not ready: executes, but hidden-history/contact-time recovery fails | Other chat |
| 3 | Not ready: strong-growth recovery improves; weak-growth identification and broader validation remain | Other chat |
| 4 | Not ready: numerical execution repaired; split-time error and calibration gaps remain | Other chat |
| 5 | Not ready: three completed cases have no agreement-gate passes; simulation memory failures remain | Other chat |
| 6 | Not ready: eight development configurations fail recovery | Other chat |
| 7 | Not ready: masked production path passes software smoke; matched BGS biological validation is missing | Other chat |
| 8 | Not ready: detection/localization and full-chromosome null calibration remain | Other chat |
| 9 | Submitted and scored: relation 0.834; meiosis count 0.842; submission 9785353 | Package: other chat; final run/UI: coordinator |
| 10 | Development passed: fixed pair-only integer-edge decoder517/585 (88.38%), macro85.29%, min50%; fresh qualification, confidence validation and229-person runtime remain | This chat |
| 11 | Not ready: expanded 64-cohort fit failed; corrected fully crossed 216-cohort panel is running as generation job10856905, fit10856907 and audit10856920 | This chat |
| 12 | Not ready: histogram repair passes software checks but fails fresh accuracy criteria | Other chat |
| 13 | Not ready: original full-width development fails; unchanged three-start full-batch refinement running as10857333 | This chat |
| 14 | Not ready: all six original confirmations fail target stability; both subsequent full-window development cases hit their fixed iteration cap; comparisons completing | This chat |

The other chat reports PR #6 merged after CI passed, covering the Task 7/12
repair snapshot with 1,022 tests on three platforms. Its Task 1 optimizer
snapshot is separate. These software checks do not close biological gaps.

## Independently verified Task 10 diagnostics

The coordinator verified the two new development cohorts against their raw
truth and prediction files, checking all manifest hashes before and after
scoring. Every method retains all 148 within-family pairs and all 14 labels.
The independent scorer's 18 focused tests pass, including missing-pair
denominators, extra-pair rejection, cohort bootstrap and calibration handling.

| Method | Correct / 148 | Accuracy | Balanced accuracy | Weighted Jaccard |
|---|---:|---:|---:|---:|
| Frozen raw-haplotype baseline | 52 | 35.14% | 42.14% | 0.2416 |
| Corrected diploid-state features | 60 | 40.54% | 46.43% | 0.3103 |
| Exact simulated diploid IBD diagnostic | 62 | 41.89% | 47.77% | 0.3215 |
| Whole-family Bonsai | 59 | 39.86% | 47.50% | 0.2982 |
| Neighbor-joining, 5 cM | 48 | 32.43% | 42.86% | 0.2207 |
| Neighbor-joining, 10 cM | 88 | 59.46% | 65.71% | 0.4288 |

Every method has at least one accepted label with zero recall. None passes
the declared advancement criteria. The two different family templates do not
support an all-label cohort bootstrap interval: some resamples omit labels,
so the scorer refuses to condition the interval on successful resamples.
Moreover, an audit found overlapping chromosome ancestry RNG seeds in the
original development generator. These results are diagnostics, not independent
confirmation. Independent *rescoring* does not establish simulation independence.

Evidence is in
`experiments/relatedness_qualification_2026_11_20/task10_development/`.
The reproducible audit is `verify_task10_development.py` in the parent
qualification directory. Original inputs remain unchanged.

The subsequent candidate used chromosome-wise sharing across phased homologs, with
explicit IBD2 handling and a matched phase-free control. Fresh training and
calibration preceded the detected-data development comparison; seed derivation
changed to separated, recorded streams. Its results below do not establish
readiness. The primary
[PONDEROSA paper](https://doi.org/10.1093/genetics/iyaf094) supports haplotype
sharing as an informative feature for close relatives, but does not establish
accuracy for this challenge's complete label set. The
[Bonsai paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8595950/) also explains
why overlapping IBD among several relatives can distinguish ancestral
lineages that pairwise likelihoods cannot distinguish. Its exact phased
diagnostic improved to 828/1,400 correct, versus 654/1,400 for the phase-free
control; these agent-reported results remain diagnostic. Same-meiosis cousin
classes remain weak. A subsequent multiway context experiment is preregistered
with disjoint training and diagnostic branch-depth multisets and a matched
pair-only ablation; it is not yet evidence of improved detected-data accuracy.

On four fresh detected development families, the phased candidate failed:
87/296 correct (29.39%), versus 81/296 (27.36%) for its phase-free control and
124/296 (41.89%) for the prior baseline. Weighted Jaccard scores are 0.2149,
0.1830 and 0.3152, respectively. The coordinator independently verified every
raw row, complete pair coverage and all manifest hashes. Raw, uncalibrated
ECE is 0.2861 for the phased model; the separately selected development
temperature is 4, the fixed grid's upper endpoint. No confirmation was run.
Paired exact-versus-detected diagnostics on identical genotypes and parentage
found substantial shifts in short-segment and homolog-balance features. The
exact-only multiway experiment remains unlaunched. A new detected-data panel is
frozen instead: 32 training, eight calibration and eight development families,
with 48 distinct branch-depth multisets and 2,160/623/585 pairs. All 14 labels
meet predeclared support minimums; reserved confirmation topologies are excluded.
The coordinator independently checked the counts, topology separation and
3,216 recorded stage seeds. A successful pilot gates dispatch of the rest,
with eight concurrent four-CPU family jobs at most. The estimated cost is
30–40 CPU-hours; the hard family-plus-fit limit is 146 CPU-hours. All 88 files
bound by its immutable source manifest independently match their recorded
hashes. This panel varies depth multisets within one outbred star-family
generator; it does not qualify arbitrary nested pedigrees, 229-person scale,
missingness or diverse ancestry.
The pilot completed in 8m51s and released the remaining 47 families. A separate
bounded exact-segment runtime preflight covers 32/64/128/229 people, without
claiming biological qualification. Its original implementation timed out at
128 people (300 seconds) and 229 people (420 seconds). A separately frozen
optimization has passed all-feature/group equivalence on eight declared
training cohorts and 26 focused tests, including exhaustive three-person
diploid graphs. At 128 people it completes in 75.7 seconds, and at 229 people in 257.64
seconds with about 1.63 GiB memory. The exact fitted-forest probability check
remains pending the model fit. Runtime success does not establish accuracy.
The proposed nested-pedigree qualification design has also passed an independent
parentage audit: all 53,764 pairs across eight small and two 229-person families
have the declared closest common ancestral couple and relationship, with all
14 labels in every family. This verifies simulation structure only; genetic
qualification remains conditional on passing development and runtime gates.
The frozen genotype path explicitly rejects
missing calls; upstream [hap-IBD](https://github.com/browning-lab/hap-ibd) also
requires complete phased genotypes. Rejection tests cannot establish robustness
to missingness, and no silent allele filling is permitted by this candidate.
The independent reports are in `task10_phase_development/` under the
qualification directory.

## Task 11 development and prospective confirmation

The frozen v3 refit reports 680/735 correct calls on four calibration cohorts
(92.52% accuracy, 85.36% balanced accuracy, 73.68% minimum label recall,
ECE 0.02768).
These are calibration results, not a release certificate. Training-only
clustering comparisons recovered 125/128 exact families with the existing
spectral method and 128/128 with either single or average linkage. Average
linkage was selected before confirmation, with the pair classifier unchanged.
Sibling pooling did not improve calibration and remains disabled.

The v3 protocol included eight untouched confirmation cohorts, four with
additional mixed-generation branches. The source archive, protocol, model,
input and result hashes are checked. The coordinator reviewed the complete
pair denominators and the separation of fitting, calibration and confirmation.
The agent reports 34 focused tests and Ruff passing in its isolated environment.

The completed confirmation **failed** the frozen criteria: 1,223/1,443 correct
(84.75%), balanced accuracy 73.82% against a 75% gate, once-removed recall
44.83% against 60%, and bootstrap lower bound 68.40% against 70%.
All 64 families were recovered exactly, with every true pair covered and no
extras. ECE was 0.038. The 160 ID/order permutations check family partitions,
not 160 independent biological cohorts or full classifier reruns. A separate
truth-free adapter smoke checks end-to-end export on a synthetic cohort.
The coordinator independently rebuilt these scores from all raw rows, checked
all 64 partitions and pair denominators, verified before/after evidence hashes,
and checked all 25 remote files bound to the candidate. The independently
computed weighted Jaccard score is 0.7701. Reports are retained in
`experiments/relatedness_qualification_2026_11_20/task11_confirmation_v3/`.
The adapter's synthetic smoke also passed: 183 pairs, no label changes and
maximum probability difference 1.67e-16. The eight exposed confirmation cohorts
are now development evidence;
they cannot serve as fresh confirmation for a revised model.

The v4 phase extension was source-frozen and independently verified against
all 231 source hashes. Its matched development ablation confirmed that the
original feature matrices stayed unchanged. On the eight exposed cohorts,
it improved to 1,254/1,443 correct and 76.85% balanced accuracy, but
once-removed recall remained 50% and half-avuncular recall 53.85%. Calibration
also regressed to 666/735 correct. These development results do
not support spending the fresh confirmation seeds. One fixed variant using
only phase segments at least 10 cM long was tested, motivated by the
paired short-segment transfer diagnostics; the acceptance criteria remain fixed.
The v5 source archive independently passes all 232 source bindings and syntax
checks. Its only changed development-script functions define the protocol,
validate its version and apply the fixed phase-length filter. IBD2 masks still
use all detected segments. Development replay binds original detector outputs
and recomputes only the 18 phase columns; it makes no claim of rerunning the
detector. Fresh confirmation, if justified, must use the complete pipeline.
The v5 development result also failed: 1,254/1,443 correct, 76.76% balanced
accuracy, once-removed recall 66/116 (56.90%) and half-avuncular recall 43/78
(55.13%). ECE is 0.0213. Both variants now pass independent raw-row rescoring,
all 64 family partitions and all 73 remotely bound input/model/result hashes.
Their weighted Jaccard scores are 0.80185 and 0.80105; both still fail the
minimum-label gate. No fresh confirmation is launched. A separate diagnostic
examines information from predicted siblings and third relatives; the
frozen phase variants remain unchanged. The primary
[DRUID study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6035284/) motivates using
siblings' combined IBD sharing, but its human results do not validate these
dog relationship labels in the presence of inbreeding.

The fixed family-context diagnostic improves the exposed eight cohorts to
1,312/1,443 correct (90.92%), 82.91% balanced accuracy, weighted Jaccard 0.85669,
ECE 0.02324 and bootstrap lower bound 77.27%. Half-avuncular recall remains
39/78 (50%), so it fails the unchanged 60% minimum. Calibration also fails
the per-label criterion. These scores were independently rebuilt from the
raw rows and archive-bound outputs in `task11_context_development/`.
Predicted sibling grouping is perfect in this diagnostic, pointing to a
classification limitation rather than missed sibling groups.

A prospectively frozen expansion adds 48 independent training cohorts to
the original 16. The combined 64 contain 32 simple and 32 mixed-generation
cohorts; each of four training folds contains eight of each. Both new topology
strata cover all nine population-size/missingness combinations. The fixed v5
pipeline and the same 36 context features, classifier and acceptance gates
are retained, with a matched base-feature ablation. All 48 new seeds are
disjoint from exposed and reserved cohorts. Generation array 10855669 and
dependent fitting job 10855708 are running/pending. All additional cohorts
must complete; no outcome-based exclusion is allowed. Fresh confirmation
seeds 2026120400–2026120407 remain unused.

## Tasks 13–14: exact official schemas and likelihood pipeline

Both official templates were verified through the authenticated Synapse UI.
Only templates were accessed; no testing or final challenge datasets were read.
The retained bytes match the official size and MD5 values:

| Task | Official template | Version | Bytes | MD5 |
|---|---|---:|---:|---|
| 13 | [bottleneck.submission.yaml](https://www.synapse.org/Synapse:syn74933357) | 1 | 208 | `d4f6e70c1fbc61055e78ab5f6fe88392` |
| 14 | [split_isolation_low_coverage.submission.yaml](https://www.synapse.org/Synapse:syn77432299) | 1 | 229 | `4b36c91baa6706d829b3ce3ff31e83e3` |

Task 13's targets are `generations` and `post_decline_fraction`; Task 14's are
`pop1_population_size`, `pop2_population_size` and `generations_since_split`.
Files and retrieval provenance are retained in
`experiments/lowcoverage_2026_11_20/official_templates/`.

The new pipeline retains genotype uncertainty through ANGSD SAF and realSFS,
then fits globally folded spectra. It avoids complete-genotype filtering and
retains the analyzed monomorphic-site exposure for population-size scaling.
Operational smokes reached realSFS, whose convergence still needs verification.
The next version adds BWA remapping of generated reads, independent Powell
checks of the demographic optimizer, explicit finite-difference settings, and
absolute likelihood-change checks. Full-length, independent-seed validation is
still required. The agent disclosed accidental exposure to public testing
answers while reading the wiki; no values are used for fitting or tuning, and
the work must not be described as wholly blinded to those public answers.

GL2 recovered much of the rare-allele mass missing under GL1, but the 1 Mb
development demographic fit still failed. A genotype-truth control, with
monomorphic exposure explicitly restored, also has a shallow likelihood ridge
and fails optimizer agreement. The full 10 Mb genotype-truth controls passed
independent arithmetic, scaling and optimizer-record checks. Task 14's relative
errors were -3.10% for split age, -9.66% for population 1 size and -3.03% for
population 2 size. Task 13's errors were -31.09% for age and +29.44% for
bottleneck fraction; age still fails the 30% criterion. Both optimizers agreed
on each full-width truth control. These controls establish neither
read-likelihood recovery nor full readiness. All 242 files bound by the
checkpoint manifest were independently verified; the manifest itself is also
bound by the archive hash. The audit is `lowcoverage_checkpoint_v1_audit.json`
in the qualification directory.

The original realSFS run missed its requested stopping rule after
5,000 iterations. Its later continuation reaches the probability-update stop,
but the full-frequency change remains above its registered stability threshold
and population 2 size is about 54% too low. A separately registered winsfs 0.7.0 development comparison
uses the same SAF evidence, independent seeds and continuation checks. Its
different stopping rule is explicitly separate from the original realSFS
acceptance criteria. The 1 Mb bottleneck comparison still fails demographic
recovery and optimizer agreement; it does not establish readiness.
The matched 1 Mb split-isolation comparison passes optimizer agreement but
still has about -30.44% population 2 error, 3.65–5.00% alternate-seed target
changes and two initial runs hitting the 20-epoch cap. A nearly unchanged
continuation alone cannot establish convergence. Before any full-width BAM
estimates, a separate development policy fixed a 50-epoch cap, the same
stopping tolerance, two seeds and continuation, less than 1% target changes,
full-frequency L1 below .001 and segregating relative L1 below .05. Recovery
and optimizer gates remain unchanged. This prospective window-estimator policy
does not retrospectively pass the failed realSFS gates.
All 16 files in the full-width source snapshot independently match their
manifest. The launch wrapper checks source, protocol and BAM/reference/popfile
audits before and after execution and binds the actual likelihood inputs.
The completed 1 Mb GL2 realSFS/winsfs comparison uses identical hashes for all
six SAF files and the same 999,993 sites. realSFS still hits its iteration cap;
the different spectra are not caused by different input files or exposure.
Full-width development array 10854727 completed. Task 13 fails independent
optimizer agreement, target stability and bottleneck-fraction recovery:
seed1 age is 2,003.70 generations, but fraction 0.28438 is 42.19% high;
alternate-seed age changes by 17.38%. Its confirmation seeds remain untouched.
A separate bounded full-batch realSFS refinement from each of the three
window-estimator spectra is registered for after Task 14, retaining original
realSFS stopping/stability criteria and all recorded failures.

Task 14 passes all full-width development gates. Seed1 split age is 4,832.37
generations (-3.35%), population 1 size 4,474.43 (-10.51%) and population 2
size 14,989.28 (-0.0715%). All runs stop before the epoch cap, both optimizers
agree, and alternate-seed target changes are below 0.785%. The coordinator
independently checked all 340 checkpoint-bound files plus the manifest, raw
spectra, actual stopping logs, optimizer records and target/gate arithmetic.
The audit is `lowcoverage_checkpoint_v2_audit.json`.

The complete Task 14 method was frozen at 19:04:18 UTC before dispatch of six
prespecified confirmation cohorts at 19:04:55. The independent freeze audit
matches all 3,194 bound remote files, including 3,166 runtime files, with no
mismatches. Original six seeds, three histories, Q30/Q20 profiles and depth
variation are retained. Preparation array 10855819 precedes inference array
10855820, each limited to two concurrent four-CPU jobs. All six cohorts and
every numerical/recovery/provenance gate must pass; one development pass
alone does not establish readiness. No challenge data is accessed.

## Published scoring formulas

The [official scoring page](https://www.synapse.org/Synapse:syn74803798/wiki/640739)
defines relatedness as a truth-class-size-weighted Jaccard index. The separate
`scripts/score_relatedness_published_metric.py` reproduces that formula from raw
rows, with missing calls retained as false negatives and extra calls as false
positives. Four hand-calculated tests and Ruff pass. This is an independent
implementation of the published formula, not verified parity with the scoring
server. It supplements the frozen acceptance criteria without changing them.

The demographic formula published on that page is the square root of the
**sum**, without division by parameter count, of squared relative errors.
Applicability to the later-added low-coverage tasks is inferred from the
general demographic-history wording; the page does not name them separately.
Prior mean-based diagnostic RMSE values must not be relabeled as that score.
Retrieval provenance is in `official_scoring_provenance.json` under the
qualification directory; complete answer/details blocks were removed before
inspection.

## Local file recovery and reproducibility

iCloud evicted some workspace sources and tests, causing reads to stall.
This chat restored only its own newly created files, either from a hash-matched
immutable compute snapshot or the exact original file-change records in this
chat. Recovery hashes are retained in
`experiments/relatedness_qualification_2026_11_20/local_file_recovery.json`.
Validation now uses isolated temporary copies and bytecode caches; shared
sources, configurations and other sessions' processes were not altered.

## Expanded Task 11 fit — independent failure audit

All 64 training cohorts completed, but the fixed context refit scored 1,231/1,443
(85.31%) on exposed development, with balanced accuracy 71.86%, half-avuncular
recall 24.36%, and cohort-bootstrap lower bound 62.65%. Weighted Jaccard was
0.7809, ECE 0.03909. The matched pair model scored 1,166/1,443, balanced
accuracy 65.13%, and weighted Jaccard 0.7296. All 64 development families
and 32 calibration families remained exact, with no missing or extra pairs.
The coordinator independently recomputed both models on calibration and
exposed data and verified the compact archive and included result hashes.
Sparse stress and fresh confirmation were correctly left unrun.
The generator is now under review for a possible coupling between ancestral
mating kind and background population size; that explanation is not yet a
causal finding. Reports are in `task11_expanded_development/` under the
coordinator qualification directory.

## Development update at 19:51 UTC

Task 10's frozen 48-family panel finished with all 48 jobs successful. Root
independently verified all 585 development pairs against eight original truth
files and sample rosters, every bound local file before/after scoring, and 24
remote truth/roster/baseline bindings. The context model achieved386/585
(65.98%, macro59.83%, weighted Jaccard0.567916), versus337/585 for the matched
pair-only model and269/585 for the prior baseline. Minimum context label recall
is zero. A supplementary NJ10 decoder reached405/585 with context and409/585
with pair-only scoring. A further fixed posterior-distance diagnostic reached
389/585 and437/585, respectively; the better pair-only version has macro67.79%,
weighted Jaccard0.633584 and zero minimum recall. These are reused development
results and authorize no qualification or submission. The cohort bootstrap
refuses intervals because580/2000 resamples lack some accepted labels. Evidence
is retained under `task10_panel_v2_development/` and
`task10_forest_distance_development/` in the root qualification directory.

Task 11's expanded fit failure is independently verified against the compact
raw archive and254 original remote bindings. A metadata audit found that
ancestry schedules were confounded with label and effective population size;
the observed failure is preserved, and a causal explanation is not established.
A new216-cohort design crosses three population sizes, three missingness levels,
three ancestry shifts and two topology settings:108 training,54 calibration and
54 development. All labels and all cohorts remain in the evaluation. Root
verified the frozen bundle, complete crosses, balanced OOF folds, source changes
and232 distinct master seeds, without claiming exhaustive uniqueness of derived
32-bit seeds. The fixed generation/fit/audit pipeline launched after that review.
Future stress and confirmation cohorts are unlaunched and remain gated.

All six original Task 14 confirmations are independently audited from717
manifest-bound artifacts, including raw spectra, stopping logs, all optimizer
starts, physical-unit conversion, provenance and complete parameter recovery.
Every cohort fails the unchanged1% target-stability requirement; one also fails
frequency L1. Median primary absolute relative error is5.354%, which does not
cancel numerical instability. The subsequent prospective full-window estimator
uses500 blocks/window500, tolerance1e-6 and cap200, keeping all other numeric
gates fixed. Both initial development runs reached the cap with log-likelihood
changes above tolerance. No fresh963xxx confirmation seeds may be consumed
under this failed policy. Complete comparisons are being retained before any
new development proposal. Root has downloaded the current official BAM archive
and reference for future final execution only; archive MD5 matches Synapse,
and a locally reconstructed FASTA index matches the official index MD5 exactly.
No Task14 final inference has run.

## Task 10 advancement at 20:01 UTC

A subsequently registered fixed-topology integer-edge diagnostic passed the
six unchanged development advancement gates for the pair-only variant. Root
independently reproduced517/585 exact (88.376%), macro85.289%, minimum label
recall50%, weighted Jaccard0.826141, complete pair coverage, baseline gain and
no cohort regression. Nine frozen remote input bindings also match. Context
reached513/585 but regressed on one family and therefore failed. The selected
pair-only rule changes only positive integer edge lengths on the inferred NJ
topology, using frozen model probabilities and at most20 improving moves/200
feasibility checks. Truth is absent from inference. All prior failures remain.

This is reused development, not fresh qualification. The agent is preparing a
new prospective source/qualification freeze, a bounded229-person runtime check,
and a separately defined global correctness confidence using calibration data
only. Any scalar confidence must be described as having no individual uncertainty
resolution, not as a root posterior; the original ECE requirement remains.
Root ambiguity and label/cohort diagnostics remain visible. No fresh nested or
229-person qualification genotype cohort has yet run.

Task13's unchanged registered full-batch diagnostic started as10857333 after
all six original Task14 confirmation jobs completed. Three4CPU tasks run at
most two concurrently,45 minutes each,9CPU-hours maximum. It uses the retained
synthetic SAFs only and does not authorize confirmation or final inference.

## Submission decision update at 20:29 UTC

The user explicitly requested more submissions if candidates are "somewhat ok".
This permits a separate practical submission decision based on credible measured
performance, while every failed strict qualification gate remains reported as a
failure. Source integrity, complete outputs, correct formatting, and truthful
provenance remain required. Official challenge rules were rechecked: each user
or team gets one final submission per task, so submitted candidates cannot be
replaced after observing their final scores.

The refreshed Task 9 dashboard has three entries. Our relationship score0.834
and meiosis score0.842 are each first, ahead of the next best0.828 and0.840.
This is the current Task 9 ranking, not an overall or closed challenge ranking.
The original Task 9 Telegram score delivery remains verified; no duplicate sent.

Task 10's selected structural pair candidate passed the229-person/26,106-pair
runtime check in16.27 seconds with complete output and verified source parity.
The fixed calibration produced483/623 correct and global confidence0.7744;
root independently verified33 local files,18 remote bindings, raw truth/rosters,
all predictions, and the constant confidence rule. This estimate has zero
individual uncertainty resolution. Eight fresh small nested cohorts are running
under the independently reviewed v3 protocol. The official final VCF downloaded
with the exact published MD5 and is being staged without inference. An immutable
final adapter is being prepared while the checks run.

Task 11's216-cohort corrected development pipeline continues. The agent is also
preparing the existing context-v1 model as a practical fallback, with90.92%
exposed-development accuracy and weighted Jaccard0.856695. Its half-avuncular
recall is50%; it has not passed strict qualification. The final choice will use
synthetic evidence and operational checks, before the single final submission.

Task 14's cap200 development independently fails only stopping-before-cap; all
other registered gates pass. Three attempts to schedule the strict cap500 plan
failed closed before scientific execution, and its full fresh-confirmation path
is closed for the deadline. Under the user's new practical submission objective,
root separately authorized two bounded cap500 development jobs on the existing
SAFs (4CPUs each,55 minutes each, latest start20:45), retaining unchanged numerical
checks and all failures. No fresh963xxx seed generation is authorized by that
practical path. A useful fixed estimate may be considered separately without
claiming the strict six-cohort confirmation passed. Task13's three refinements
all failed stopping, recovery, and continuation stability and remain unsuitable.


## 20:52 UTC — practical submission preparation and independent checks

Task 9 remains first of three posted entries on relationship (0.834; next 0.828) and meiosis count (0.842; next 0.840). Its score message was delivered to Telegram; do not duplicate it. The public Task 10 table currently has one competitor at 0.997 for both metrics. There is no overall challenge ranking inferred from these task tables.

Task 10 final input has 263 people and requires all 34,453 pairs. Adapter v2 source bindings and eight-cohort map/preparation identity checks passed. The 25% whole-site thinning diagnostic preserved every one of 585 predictions. Fresh small-panel strict qualification FAILED: one initial NJ tree had no accepted root, leaving 120/1,552 predictions missing. Root independently verified all raw rosters/truth/predictions, 76 local manifest files and 98 remote hashes. Structural output scored 1,432/1,552 counting missing as wrong; pair-only scored 965/1,552 and baseline 730/1,552. Seven successfully inferred families were completely correct. No original large run was launched. A separate explicit pragmatic structural-if-valid/pair-only-on-exact-no-root policy is being frozen and tested before two untouched large cohorts; old strict failure is preserved, and observed small families are now development evidence.

Task 11 first final run 10858497 failed before classification because detected segments exceeded frozen map endpoints. Source-reviewed deterministic inclusive map-domain staging passed a synthetic spike/identity check and 203-pair full-model parity; root checked 43 remote bindings. Retry 10858887 uses the same frozen context candidate and distinct output directory, preserving the failed run. Final choice still awaits the matched 54-cohort comparison of old versus corrected models.

Task 14 separate pragmatic cap500 development is running. Q30 is complete and passes all original case-level gates; Q20 is pending. Both complete cases and root raw audit are required before final inference. The separate practical policy was frozen before outcomes, retaining fixed seed1 primary and 5% target sensitivity, with strict failures still recorded. No final Task14 inference yet. Root verified all 941 checkpoint-v5 files and reconstructed all six Task13 refinement estimates/spectra/stopping checks: decline-fraction error exceeds 50%, so Task13 is closed without submission.

The other active chat retains ownership of Tasks 1–8 and 12; latest read-only snapshot reports Task3/4 PR9 and Task2 export fixes under review. No duplicate work has been started in those scopes.


## 2026-11-20T21:07:51.552100+00:00 — user explicitly prioritizes all final entries by deadline

The user now requires submissions for every remaining task regardless of validation status, with Tasks11and14 first. Validation failures must remain disclosed but are no longer submission vetoes. Root retains actual-input provenance, output format, correct queue and one-final-entry checks; no fabricated values or final truth. The user explicitly authorized messaging the other active chat; deadline-first instructions were delivered for its Tasks1–8and12, with candidate handoff requested by22:30UTC and central uploads owned byroot.

Task10 final job10859028 completed0:0 in9m17s; two prospective229family checks continue. Task11 existing-context output audited265pairs/68people while all54matched model selection completes. Task14 first final10859086 failed beforeestimation on official popfileheader; exactheaderremoved in separateinputstage, all40BAM/index files unchanged. Retry10859136 started21:01:56UTC. Task13 reopened under userpriority: best existing continuedwinsfs selected by prior whole-target error and stopping/optimizerreliability, still41.21%declinefractionerror; finalarchive/reference downloadedandverified, remote staging running.


### 21:27UTC deadline execution update

Task10 submitted successfully as individual to GHIST2026Relatedness2Final; uploaded syn78780148v1 SHAee9bace7bfc086281be0b243c5b440db9542277466b4722a14cb32523ef059be. UI confirmation retained; score/receipt ID pending dashboard. Root independent packaging audit263people34453pairs passed; final structural run had no fallback. Task11 selected correctedcontext after fullmatched54 root rescoring (.835006230874145WJ vs oldcontext.7482403619144973) and933remoteclosurehashes+129pair alias/order parity; strictper-label stillfailed. Finaljob10859650 launched21:24UTC. Task13 fixedcontinued finaljob10859453 launched21:15UTC; normalized exact organizerheader only, original40BAM/BAIunchanged. Task14 retry10859136 stillrunning. Otherchat Tasks1,3,6,12 ready; rootindependentaudit underway, upload remainscentralized. Latestuserdeadlinepriority overrides strictscientificsubmissionveto, notformat/provenance requirements.


## Scored closeout, 21:44 UTC

Tasks 11 and 14 are SCORED (submissions 9785361 and 9785360). Exact scores, receipt and candidate hashes, evidence scope, and preserved scientific limitations are recorded in [the scored closeout](continuation_scored_closeout_2026_11_20.md). Task 10's first submission was INVALID and is being corrected; the coordinator continues remaining submissions.


## Final scored closeout, 22:01 UTC

All 14 tasks now have valid server scores. [The consolidated scored summary](final_scored_summary_2026_11_20.md) lists every valid submission and score; its linked machine audit verifies exact candidates and email evidence. All current receipts are SCORED, and six initial invalid attempts remain preserved separately. Scientific qualification records remain unchanged. No further submission is required for these 14 tasks.
