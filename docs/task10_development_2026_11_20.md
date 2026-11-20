# Task 10 development toward complete-family inference

This work changes only new development files and isolated experiments. It does
not change production models, shared converters, existing validation artifacts,
or work assigned to other chats. No challenge answers, challenge testing data,
final data, or submissions enter these experiments. Task 10 remains unqualified
for submission.

## Measured failures and bounded repairs

On the original 18-person family, the frozen classifier scores 56/153 (36.60%)
with exact ancestral diploid IBD segments and 58/153 (37.91%) with detected raw
haplotype segments. Thus the principal failure cannot be attributed to the
IBD detector or repaired Bonsai adapter. Many wrong labels have the same total
meiosis count as the correct relationship.

A new development adapter turns detected haplotype edges into diploid IBD
states, requiring a genuine two-haplotype matching for IBD2. This matches the
segmentation in the frozen classifier's training input. It scores 63/153
(41.18%), including all nine sibling pairs instead of five. All pairs come from
the genotype roster, including zero-sharing pairs.

A whole-family candidate removes threshold branch partitioning: Task 10's
single family is reconstructed together under two fixed Bonsai settings,
maximum add degree/relationship types 5/2 and 5/3. Corrected haplotype conversion
uses 10 cM segments; consensus retains the previous replicated-topology and
flat-model meiosis guard. This gives 70/153 (45.75%).

A separate neighbor-joining candidate fits an unrooted tree to kinship-derived
meiosis distances, rounds edges, and enumerates valid roots. Roots are ranked
using a composite sum of frozen pair-model log scores, not an independent
likelihood or calibrated posterior. On the original family it recovers all 153
meiosis distances but only 81/153 (52.94%) exact labels. Four admissible roots
remain; the leading two scores differ by only 0.0144. No midpoint, age, generation,
or truth-selected root prior was introduced. Current diagnostics separately
report root-invariant labels and ambiguous label sets; their uncertainty is
explicitly uncalibrated.

Two fresh mixed-depth diagnostic families cover all fourteen accepted labels
across 148 complete pairs. One has sixteen people in sibling branches at root
depths 2,2,3,3,4,4,5,5; the other has eight at depths 1,2,3,4. No unsupported pairs
are discarded. Seeds 2026124010 and 2026125010 were fixed before running.

| Method | Correct pairs / 148 |
|---|---:|
| Frozen raw haplotype baseline |52|
| Diploid-state adapter |60|
| Exact ancestral IBD oracle |62|
| Whole-family Bonsai consensus |59|
| Neighbor joining, 5 cM |48|
| Neighbor joining, 10 cM |88|

**Design limitation discovered afterward:** the reused founder generator added
100 times chromosome number to the cohort seed. Because these two cohort seeds
differ by 1,000, some chromosome random seeds coincide. Different chromosome
lengths/maps do not establish independence of those streams. These results are
retained diagnostic evidence and must not be called independent confirmation.
An independent rescore refers to the separate verification procedure only.

The [raw scorer inputs and SHA manifests](../experiments/task10_development_2026_11_20/scorer_inputs/)
were independently checked against retained truth/prediction rows. Every variant
includes 148 pairs and all 14 labels. None passes the prospective advancement
gates: complete coverage, 80% exact, 70% macro recall, 50% minimum-label recall,
10 percentage-point improvement over baseline, and no cohort regression. These
are development utility targets, not an official challenge cutoff or sufficient
scientific release criteria. The [independent reports](../experiments/relatedness_qualification_2026_11_20/task10_development/)
also identify sparse class coverage across only two cohorts.

## Prospectively fixed phased candidate

The next bounded candidate retains homolog information discarded by the flat
model. It uses raw and diploid-state segment summaries plus chromosome-wise
minor-homolog sharing fractions and lengths, with IBD2 regions masked. Duplicate
and overlapping edges are unioned. Features are invariant to sample-order swaps
and independent chromosome homolog relabeling. A control model retains the
same explicit IBD2 summaries and omits only minor-homolog features.

This design is motivated by the haplotype-score method in
[PONDEROSA](https://pmc.ncbi.nlm.nih.gov/articles/PMC12341879/), whose published
validation does not establish performance for this full 14-label task. The present
extension therefore requires its own independent evidence; it is not a change
to a PONDEROSA package dependency.

The [fixed protocol](../experiments/task10_development_2026_11_20/phase_protocol.json)
predeclares 14,000 new exact phased Ped-sim pairs for fitting, 1,400 additional
exact pairs for a diagnostic, and a fixed ExtraTrees classifier with 500 trees,
minimum leaf size 5, balanced class weights, and all features available at each
split. There is no hyperparameter search. The otherwise identical control
excludes the homolog-sharing features.

Four new detected-genotype development families are separately reserved for
one scalar temperature calibration. All three previously examined families are
excluded from fitting and calibration. Confirmation seeds are reserved and
remain unexecuted until development gates pass and source, model, and scalar
calibration are frozen. Complete simulated genotypes remain phased, error-free,
and missingness-free; founders use a simple constant population size 10,000,
mutation rate 1.25e-8, chromosome-average recombination and common-marker
ascertainment. These do not represent the entire challenge distribution.

For the new batch, ancestry, mutation, meioses, and hypothetical dropout seeds
are derived by NumPy SeedSequence from cohort, chromosome and stage identifiers.
All derived seeds are recorded before generation. Tests verify determinism and
absence of collisions over all planned cohorts and stages. This repairs future
simulations; the historical diagnostic records remain unchanged.

Completed work: 10851850 (fixed phased-model training), 10851851 (four new genotype
families), 10851870 (phase/control predictions), and 10851979 (fixed calibration). All source snapshots and large synthetic artifacts live under Betty
`experiments/task10-development-20261120`; no running or completed snapshot is
modified by later local edits. Original diagnostic jobs 10851292, 10851375,
10851374, 10851526 and 10851595 are complete. Job 10851492 failed before inference
because the remote environment lacks pytest, and was replaced by 10851526.

Local unit checks are run against exact source copies outside the iCloud-backed
workspace when shared imports become dataless. They cover complete genotype
pair enumeration, phased interval matching, all 14 pedigree labels, genotype
selection, immutable-input checks, additive-tree recovery, phased invariance,
IBD2 masking and random-stream separation. No test count or software execution
success establishes scientific submission readiness.

## Phased candidate outcomes and next diagnosis

The 1,400-pair exact phased diagnostic gives 828 correct labels (59.14%) for
the phase model and 654 (46.71%) for the phase-free ablation. The clearest
improvement is avuncular once removed, from 37/100 to 96/100. Same-meiosis
all-cousin distinctions remain poor even with exact ancestral segments.

On the four new detected development cohorts, all 296 pairs and all fourteen
labels are retained. The phase model scores 87/296 (29.39%), the phase-free
control 81/296 (27.36%), and the original frozen baseline 124/296 (41.89%).
Phase macro recall is 39.11%; minimum class recall is zero. All utility gates
except label coverage fail. Scalar calibration chooses the largest predeclared
temperature, 4.0; this cannot change the predicted label or repair poor ranking.
Avuncular once removed improves from 1/8 to 8/8 relative to the ablation, but
distant cousin performance collapses. No confirmation seeds were consumed.

[Raw phase scorer inputs](../experiments/task10_development_2026_11_20/phase_evidence/phase_calibration/)
and [SHA manifests](../experiments/task10_development_2026_11_20/phase_scorer_inputs/)
retain separate phase, phase-free, and original baseline reports. Development
probability tables are uncalibrated; the scalar selection is recorded separately.

A [bounded multiway protocol](../experiments/task10_development_2026_11_20/multiway_protocol.json)
reserves 200 training and 40 diagnostic whole families with disjoint sorted
branch-depth multisets. Its proposed features use only inferred sibling cliques
and shared homolog components across two or three groups. Duplicate edges are
merged; impossible components joining both homologs of one person are counted
as conflicts and excluded. Same-training pair-only ablation and a fixed 20%
edge-dropout diagnostic are specified. This experiment is not yet launched: the
large exact-to-detected shift is now the immediate diagnosis. A paired rerun of
the four new gene drops will retain exact homolog identities, verify identical
genotypes and parentage, and compare fixed feature-group swaps without fitting.

The paired diagnostic (10852103) completed with byte-identical decompressed
genotypes and parentage for all four families. On these same 296 pairs, the
phase model scores 145 with exact ancestral sharing versus 87 with detected
sharing; the phase-free model scores 118 versus 81. Replacing only detected
phase scores with their exact counterparts restores 115 correct predictions;
replacing only raw spectra restores 111, and only state spectra restores 98.
For the phase-free model, replacing state spectra alone restores 118.

The observed shift includes short/background segments and minor-homolog
sharing. Fourth cousins have a median zero exact 2–5 cM segments versus 12.5
detected segments; median total sharing above 5 cM rises from 11.6 to 29.0 cM.
Their minor-homolog fraction rises from zero to 0.136. Second cousins have
median 2–5 cM counts 2.5 versus 15.5, while median total above 10 cM remains
134.24 versus 134.34 cM. Extra detected sharing is not automatically a caller
false positive: coalescent founders can share older ancestry absent from the
gene-drop ancestry record. An exact-only training distribution omits this
background and detector behavior. The current forest cannot advance, and the
unlaunched multiway experiment needs a transfer design that addresses this.

[Paired diagnostic scorer inputs and manifests](../experiments/task10_development_2026_11_20/paired_phase_scorer_inputs/)
retain all feature-swap outcomes. These are diagnostic interventions on reused
development cohorts, with no independence or confirmation claim.

## Detected whole-family panel v2

The new [v2 protocol](../experiments/task10_development_2026_11_20/detected_panel_protocol_v2.json)
trains both context and pair-only forests on detected features from 32 complete
families (2,160 pairs), calibrates on eight (623 pairs), and evaluates development
on eight further families (585 pairs). Every split covers all fourteen labels.
All 48 sorted branch-depth multisets differ across splits and exclude both old
reserved confirmation topologies. Each family has its own coalescent founders.
The protocol records 3,216 stage seeds, including reserved dropout streams;
a broader Cartesian audit checks 4,416 chromosome/stage combinations for
uniqueness and separation from the eight previous phase development/confirmation
families. This genotype-only design simulates no reads.

The primary context forest and matched pair-only ablation have fixed 500-tree
settings. One scalar temperature per model is chosen only on calibration data,
then model and calibration hashes are recorded before development outcomes are
loaded. The old confirmation seeds remain untouched. No threshold search or
truth-based sibling groups, roots, or generation features are introduced.

All imported sources, binaries, maps, and baseline/sibling models are copied
into a read-only Betty snapshot, `detected-panel-source`. Its manifest hash is
`7d4c00c9cfd9d197fce722e5de566078336ff0e6c2d7a39a84504e28a3bd9489`;
the protocol hash is
`e6e5ce0a71758efea55bd33cffcd80fa0bf334a50468b7ea2d7cc5619ba58f5f`.
Every job verifies this manifest. Exact homolog and detected outputs are both
retained for later paired diagnostics. Pilot job 10852393 tests the first
training family before the remaining array starts. The observed resource
estimate is 30–40 CPU-hours; the hard cap is 144 CPU-hours for 48 families at
45 minutes and four CPUs, plus two CPU-hours for fitting. Concurrency is capped
at eight families (32 CPUs, 176 GB). These resource bounds are distinct from
scientific advancement gates. Current isolated unit checks: 21 passing.

The pilot completed in 8 minutes 51 seconds. Array 10852504 runs the remaining
47 families with concurrency eight; fit/calibration/development job 10852569
depends on all simulation jobs succeeding. Before any held-out outcomes were
observed, a separate frozen runtime preflight was dispatched as job 10854208.
It profiles the unchanged feature implementation on exact synthetic segments
for 32, 64, 128 and 229 declared people, with individual limits of 180, 240,
300 and 420 seconds. Its total allocation is at most two CPU-hours. Protocol
SHA256: `babed466fbce3ec443dabc89c97b98a0b927de831df96259e1b5003c79e1ff26`.
This measures implementation scaling; it does not qualify detector accuracy.

A prospective nested-pedigree qualification design is also frozen, but genetic
qualification has not been dispatched. Four nested branching structures contain
16, 20, 20 and 24 people; two independent genetic replicates of each, plus two
229-person nested families, yield 53,764 declared pairs. Every family covers all
fourteen labels. Numerical accuracy, calibration and runtime gates are retained
in `nested_qualification_protocol.json`, SHA256
`8455e205bbd2c257225a1c0bd07c365b78f0e6260211d0e7ba7626986c94e18a`.
All 6,072 Cartesian stage/chromosome streams across the new roster and previous
panel/development/reserved-confirmation cohorts are unique. Ped-sim dry-run job
10854846 checks only the proposed parentage. The qualification runner and model
hashes must additionally freeze before dispatch, and all current development
gates must first pass. The proposed hard allocation cap is 48 CPU-hours, not a
runtime estimate. Current unit checks: 22 passing; Ruff passes.

The supported domain remains complete, chromosome-wide phased genotypes from
outbred pedigrees with unrelated spouses. The frozen hap-IBD pipeline rejects
missing alleles. Checking rejection under 5% and 20% dropout is an input-contract
check, not evidence of biological robustness. Any future imputation or site
deletion policy requires separate detection-matched development. At 229 people,
dropping every site with any missing genotype would retain only about
`0.95^229 = 0.0000079` of sites under independent 5% genotype missingness.
Neither inbreeding, ancestry shifts, sequencing errors nor phase errors are
qualified by the proposed clean-input experiment.

Structural job 10854846 stopped before generation because this Ped-sim build
requires an explicit crossover model even with `--dry_run`. The failed snapshot
and log are retained. A separate validation snapshot adds the already frozen
interference model; job 10855084 verifies all ten proposed structures and all
53,764 pair labels against Ped-sim's actual parentage. No genetic data were
generated. `nested_structure_evidence/` retains parentage, definitions, declared
rosters and truth tables; `nested_structure_manifest.json` binds their hashes.

The original scaling preflight completed with 32-person features in 4.80 seconds
and 64-person features in 23.17 seconds, but timed out at 300 seconds for 128
people and 420 seconds for 229 people. Runtime evidence is retained separately
in `runtime_preflight_evidence/`. A bounded CPU profile confirms that 128 people
produce 115 inferred groups, with pair features/grouping taking 2.76 seconds;
three-group component sweeps dominate the remaining work. This is a material
scaling failure, independent of the pending scientific development outcome.

The separately frozen optimization candidate `fast_multiway_features.py`
preconverts intervals once and limits three-group sweeps to intervals where
at least two types of cross-group edges overlap. Outside those intervals, a
component cannot span all three groups; every conflict is already present in
a two-group induced graph. Existing pair-conflict intervals are unioned with
conflicts found during the restricted three-group sweep. The genotype policy,
thresholds, inferred groups and features are unchanged. The original frozen
candidate and running jobs are not modified.

Equivalence tests include all 4,096 three-person diploid edge graphs, twenty
randomized multiperson examples with duplicate/overlapping edges and nonlinear
maps, and full-feature sample/order/homolog invariance. The initial implementation
test caught accidental base-pair adjacency merging on genetic coordinates; the
new implementation uses an explicit continuous-coordinate union. All 26 tests
and Ruff now pass. Protocol SHA256
`431cbe63c94c32f3d10fa77e8260075d48a46c22535a19366b363462eee278f8`
predeclares feature equivalence on training families 0–7 and repeat scaling on
the same four inputs and limits. Job 10855790 has a hard two-CPU-hour allocation.
Dependent job 10855853 checks exact float32 model-input and serial forest
probability equality after model freezing, using training data only. Its hard
cap is one CPU for five minutes. Optimization success would remove a runtime
blocker; it supplies no biological qualification.

Optimization job 10855790 completed. All features and inferred groups agree for
eight predeclared training families (560 pairs), with maximum scaled numerical
error `1.02e-15`. Optimized feature times are 1.56, 10.26, 75.71 and 257.64 seconds
for 32, 64, 128 and 229 people. The 229-person calculation covers all 26,106 pairs,
infers 115 groups and uses approximately 1.63 GiB peak RSS, within the unchanged
420-second cap. These are exact ancestral segment inputs. Full detected-data
229-person runtime, model probability equivalence and biological qualification
remain outstanding; the result does not establish submission readiness.

The qualification runner and nested MCAR missingness-contract helper are prepared
as unlaunched drafts. The runner requires all development, independent rescoring,
feature equivalence, model equivalence and runtime gates before creating data.
It verifies actual simulated parentage against the declared nested labels and
retains the original baseline on identical cohorts. Two scheduler drafts preserve
the predeclared 45-minute/four-CPU small-family and 90-minute/eight-CPU large-family
caps. No independent genetic qualification has launched.

The missingness helper generates nested 5% and 20% whole-genotype masks, using
per-chromosome recorded seeds and canonical sample order. Its test verifies
determinism, threshold nesting and sample-column-order invariance. The default
preparer must reject the missing alleles without publishing an output VCF; this
is explicitly a rejection-contract check, not robustness. The current owned test
suite has 27 passing tests. The historical 26-test optimization snapshot is
retained separately with its original protocol-bound SHA256.

A live runtime-environment audit before fitting matched the frozen Python and
package versions. The interpreter environment is shared and is not a hermetic
snapshot; a second version audit is required after fitting.

Before observing v2 outcomes, a single supplementary structural diagnostic was
frozen (protocol SHA256
`2b91ebb4721a45d2301cdb9cbc42c677aa9fc573719d4cbb7fe1a65f44cce868`).
It uses the existing inferred NJ tree at the fixed 10-cM distance threshold,
ranking admissible roots using frozen detector-trained context probabilities,
with an otherwise identical pair-only ablation. There is no midpoint, age,
generation or true-topology prior. Root scores remain composite and uncalibrated;
all valid-root ambiguity is retained. Job 10856408 is dependent on fitting,
limited to two CPUs for ten minutes, and automatically skips the comparison if
the original v2 gates pass. This diagnostic does not alter those gates, consume
independent qualification data or establish readiness.

The 48-family detected panel and fit completed successfully. The family jobs
used 28.47 allocated CPU-hours (elapsed time times allocated CPUs), below the
30–40 CPU-hour planning estimate and the 144 CPU-hour hard family cap. The
independent source/version audit after fitting still matches the recorded
environment. The eight development families contain 585 pairs and all 14 labels.
Context predicts 386 correctly (65.98%), compared with 337 (57.61%) for the
matched pair-only ablation and 269 (45.98%) for the original baseline. Context
macro recall is 59.83%; minimum label recall is zero, including 0/8 fourth cousins
and only 4/36 second cousins. It regresses against the original baseline in
cohort 041 (13/45 versus 14/45). Both scalar temperatures selected on the separate
calibration families are 1.0. The original exact, macro, minimum-label and
no-cohort-regression gates fail. No reserved confirmation or nested/229-person
genotype qualification has launched. These development improvements do not
establish submission readiness.

The optimization's final model check also completed: all 560 pairs in training
families 0–7 have bit-identical float32 model inputs and serial prediction
probabilities for both frozen final models. Combined with the 257.64-second,
1.63-GiB 229-person exact-segment runtime, this clears the separate implementation
equivalence and synthetic runtime checks. It does not qualify detected-data
229-person accuracy or runtime. The manifest
`optimized_runtime_and_model_manifest.json` binds both reports.

The conditional fixed NJ10 comparison ran after the primary model failed.
Context-root scoring predicts 405/585 (69.23%; macro 66.31%) and pair-only-root
scoring 409/585 (69.91%; macro 66.78%). Both have zero minimum recall, including
0/56 second cousins twice removed. There are one to six valid roots per cohort;
none uses the no-valid-root fallback. Predictions remain exploratory and root
scores uncalibrated. A subsequent explicitly oracle diagnostic finds that even
choosing the best valid root using truth reaches only 417/585 (71.28%). Inferred
unrooted meiosis distances are correct on 431/585 (73.68%). This separates a
substantial tree-distance/topology error from the smaller root-ranking error;
the oracle labels are never exported as candidate predictions.

Independent-rescoring exports include every raw truth/prediction table, original
roster, prospective protocol, calibration/model freeze and SHA256 manifests.
Retained truth, roster and baseline bytes match both simulation-output hashes
and fitting-input bindings. Parent verification is requested for all three v2
variants and both supplementary tree variants. All earlier failed policies,
artifacts and frozen sources remain unchanged.

The next fixed diagnostic replaced the raw kinship-derived distance with the
frozen forest probability-weighted mean meiosis count. Source/protocol were
frozen before that comparison (protocol SHA256
`dbc38a3071b17883653580c7ceaf346af735430a81d5a673d256eca4ae260ea1`),
using explicitly reused development families. Job 10857019 completed in eleven
seconds. Context scored 389/585 (66.50%; macro 58.59%) and pair-only 437/585
(74.70%; macro 67.79%). Both retain zero minimum recall; no qualification is
launched. A separate oracle diagnosis places the best possible root on the
pair-only inferred trees at 461/585 (78.80%), with correct meiosis distances for
469/585 (80.17%). The corresponding context values are 397 and 405. Oracle
root choices remain diagnostic only.

One additional bounded integer-edge diagnostic is prospectively frozen after
these failures (protocol SHA256
`a5434536fcf074fba1d4f8098649c8fe0579f248abcd2407c8596ebfae35bcae`).
It retains each variant's inferred mean-distance topology, compresses only
degree-two nodes, and changes one positive integer edge length by one per move.
The objective sums log model probabilities marginalized over labels with the
current pairwise meiosis count. It accepts at most twenty improving moves and
performs at most two hundred feasible-root checks. The final full-label root
score is unchanged. Root pruning only excludes positions that necessarily
violate the accepted-label maximum of five meioses; tests compare the pruned
and unpruned valid roots. No true topology, root prior or outcome-dependent
exclusion is used. All 32 tests and Ruff pass. Job 10857310 has a hard two-CPU,
ten-minute cap. Root ambiguity remains conditional on chosen edge lengths and
does not cover alternative length assignments or topologies. The original gates
and untouched qualification data remain unchanged.

Job 10857310 completed in seventeen seconds. The pair-only integer-edge
candidate scores 517/585 (88.38%), macro recall 85.29%, and minimum label recall
50%, passing all six unchanged development advancement gates. Context scores
513/585 (87.69%), macro recall 84.18%, and minimum recall 50%, but fails the
individual-cohort regression gate (cohort 041: 13/45 versus baseline 14/45).
Each fit accepts zero to two edge moves; neither cap is reached. These are
reused development outcomes following several retained failed methods, not
independent confirmation. Raw rows, inferred trees, source bindings and
manifests are exported under `integer_edge_scorer_inputs/` for independent
verification. The old nested protocol targets the earlier context candidate;
it cannot silently qualify this newly selected pair-only structural candidate.
Qualification remains unlaunched, including the 229-person genotype runs and
confidence-calibration gate. Root composite scores remain uncalibrated.

Parent independently verified both integer-edge exports, original bindings and
all development gates; pair-only weighted Jaccard is 0.8261413. A separate
integrated pair-only adapter avoids the unused context computation. Operational
preflight 10857530 completes all 26,106 pairs of the retained 229-person exact
ancestral segment fixture in 16.268 seconds: loading 0.707, features 7.424,
forest prediction 1.282, and structural inference 6.856 seconds. Python peak
RSS is 783,072 KiB; final Slurm batch MaxRSS is 2,249,292 KiB. Both fit the fixed
10-GiB cap. All 560 retained training pairs have bit-identical float32 inputs
and serial forest probabilities; all 585 reused development structural labels
are unchanged. Parent independently verified 43 input bindings and nine runtime
source bindings. This operational fixture is not detected-data qualification.

The separate prospective v3 protocol retains the ten previously declared
nested rosters/seeds, all 53,764 pairs, and original accuracy, recall, coverage,
baseline-gain and ECE gates. Eight small nested families must pass before the
two 229-person genotype qualifications launch. It adds a fixed global
correctness calibration from all 623 pairs of the existing eight calibration
families: `(correct + 1) / (pairs + 2)` for the selected structural label, with
uniform residual probability over the thirteen alternatives. This has zero
individual resolution and is not a root posterior. Reliability by true label,
predicted label and cohort remains visible; ECE must be at most 0.10 for the
small aggregate, each large family and pooled qualification. Every root
ambiguity is retained conditional on the inferred topology and edge lengths.
Calibration policy and full source freeze precede calibration outcomes and
qualification genotypes. All 34 tests and Ruff pass. At this documentation
update, calibration and qualification remain unlaunched pending source review.

The final v3 protocol is `9b9ebaa12e95ce0cd85c0be8056c33ede2a7dc1af5c717594c21383cd46d5803`,
with remote qualification source manifest `d59da1b628ff57bc7a7bcccd41503e44408f1c5f76d62b1666b5d1fcf6864036`.
Task11 performed an independent read-only source, import-origin, model-binding
and ten-roster review before launch. Calibration job 10857997 completed with
483/623 correct; the fixed global confidence is 484/625 = 0.7744. Parent
independently verified all raw rows, 33 local bindings and 18 original remote
bindings. This scalar cannot express individual uncertainty. Eight fresh small
nested qualifications launched as array 10858090, concurrency four, after the
confidence and review bindings were made immutable. No large qualification had
launched at 20:41 UTC. Reserved original confirmation seeds remain untouched.

The user subsequently relaxed the practical submission target. A separately
named practical scaling declaration, SHA256
`3efb0c2e0333a905fd1ba2009be76990e51fde55314d7ac74e6225b1e9e8ea49`,
was recorded before observing the small nested outcomes. It permits the same
two fixed 229-person cohorts after all eight small cases and an independent
raw audit, if exact accuracy is at least 70%, macro recall at least 60%, weighted
Jaccard exceeds the matched baseline, and every output/resource check passes.
This does not change the strict qualification criteria. A separate new runner
preserves the original strict result, explicitly labels practical output, and
rejects actual starts at or after 21:50 UTC. It retains the model, seeds,
generation, detector and inference exactly. Its frozen source manifest is
`80bef98356b39c8a002ddec0e8f0e5fa7bd53868120f96551fb14ca154ba580e`;
no practical large run has been dispatched at this update.

The inference-only handoff uses the same frozen candidate and original sample
identifiers, validates all 229 samples/26,106 pairs, and atomically writes the
three-column relation CSV. Parent owns actual final-data access, execution,
packaging, selection and submission. A deterministic inclusive crop to the
frozen genetic-map intervals excludes and counts outside-domain markers without
extrapolating cM or changing calls; training simulations had no outside-domain
markers. The minimum 100 markers per chromosome is the operational caller
minimum, not evidence of equal marker density. All import origins, package
versions, source, map, model and confidence are bound. Handoff protocol SHA256
is `29ad1f1262235f170585a777cca9338bbe6a920a3a8957e417d0a0328877a2b4`;
manifest is `46fd29342db0bedb0c91b4a44721a83106b16eb669d32f7c29a92a4c3fa1be01`.
Four contract tests pass, covering exact endpoint behavior, unchanged calls,
identity crop, rejection without partial publication, complete 229-person CSV
coverage and deterministic site thinning. Scheduled job 10858557 verifies crop
and preparation byte identity after decompression on all eight retained
whole-family development VCFs; its outcome is pending at this update.

Because official phased/error-free inputs may contain missing genotypes, the
strict path still rejects missing alleles. A separate pragmatic complete-site
option was declared before final inspection, with no allele filling and at
least 75% of in-map biallelic variant markers retained on every chromosome.
Its frozen declaration SHA256 is
`3648c748c8ce1663a1a11167671e1c88890c3a235fad7c1eb68caa48213b9abc`.
The fixed diagnostic uniformly removes exactly floor(25% of sites) per
chromosome from all eight reused development VCFs, using recorded SeedSequence
streams. All jobs in array 10858495 completed in 20–27 seconds; the method is
unchanged and no fitting occurs. The fixed assessor finds 517/585 correct,
85.2891% macro recall, and all fourteen label recalls identical to intact data,
passing the declared maximum three-percentage-point accuracy/macro losses and
minimum-class preservation gates. Parent independent raw verification remains
required before this optional path is used. This MCAR whole-site thinning
check is not independent confirmation or validation of structured missingness,
allele-dependent missingness, phase error or all missingness rates. Strict
missingness qualification is unchanged. The simulator's human mutation rate
1.25e-8 also differs from the official stated 1e-8; this domain limitation is
retained explicitly in the handoff.

Parent's read-only official final preflight subsequently found **263** unique
samples, rather than the older assumed 229. No individuals may be discarded.
The minimal inference handoff v2 therefore requires `--expected-people` from
the complete root-preflight header, checks exact equality, and emits every
`n*(n-1)/2` pair: 34,453 for 263. No feature, detector, model, confidence or
structural rule changes. Both 229- and 263-person CSV coverage/rejection tests
pass (five tests total). V1 remains immutable, and the 229-person operational
and synthetic qualification evidence is not relabeled 263-person evidence.
V2 protocol SHA256 is
`5ce7ed75ae31c10ee115d26e10d03edc98f54119b914443ccfd1086ea6db793c`;
manifest SHA256 is
`5daf3871628bd87dd309b2cfb7231620ee00dec307956ee45a1cbc7f2dc0ee20`.
Parent retains actual final inputs; this agent has not accessed or modeled them.

The complete small-stage assessment is **INCOMPLETE_QUALIFICATION_FAILED**.
Cohort 000's inferred initial NJ tree has no accepted root; the frozen candidate
raises its declared error and emits no structural labels for its 120 pairs.
The other seven cohorts complete, with all 1,432 emitted structural labels
correct. Counting every missing label as wrong yields 1,432/1,552 = 92.2681%;
this is not full output coverage. The strict assessor preserves the failure
and does not gate on a selected seven-cohort subset. The separately proposed
practical stage also fails its complete-output prerequisite, so no large
qualification or practical scaling was dispatched on these results. Original
logs, raw truth/rosters/predictions and hashes are retained in
`structural_small_evidence/`; manifest SHA256 is
`c827bad87662be183460e2ce28bc51fc4cc7de2151a14f75d33d48c01ed9150f`.
The frozen pair-only alternative emits all 1,552 pairs and scores 965 correct
(62.1778%); for the failed structural cohort it scores 78/120. A mechanical
structural-if-valid, pair-only-on-no-root policy would score 1,510/1,552
(97.2938%) on these now-observed cohorts, but that is a post-hoc diagnostic of
a different policy, not a qualification result or implemented fallback.
The root was informed immediately; any new policy requires its own explicit
freeze and prospective validation while retaining the original failure.
Comparator-only job 10858815 computes the missing cohort-000 baseline in a
new directory and exports all eight cohort scorer inputs, leaving structural
cohort 000 empty rather than fabricating predictions.

Crop/preparation parity job 10858557 passes all eight retained development
VCFs. Every original-versus-cropped decompressed hash and every old-versus-new
prepared decompressed hash matches, with zero excluded markers. This verifies
that the declared map crop has identity behavior on the training domain; it is
not new biological evidence. V2 only changes roster-count validation, so its
crop and genotype preparation are the same verified functions.

Parent independently verified the incomplete small-stage evidence, including
76 local manifest files, 98 original remote hashes and every one of the 1,552
truth/roster pairs. The completed baseline is 730/1,552 (47.0361%), pair-only
965/1,552 (62.1778%), and structural 1,432 correct with 120 missing. Both old
strict and old practical gates remain false. The user-authorized pragmatic
follow-up is a new explicitly declared conditional policy: retain the frozen
structural result if inference succeeds; only the exact initial-NJ-no-accepted-
root ValueError permits the frozen pair-only output for the entire cohort.
Every other exception propagates. It uses no truth and no per-pair selective
replacement. Complete ordered pair coverage, roster and input-row digest are
recorded along with selected method, fallback flag and exact reason.

This new policy is `structural_fallback_policy.py`, frozen under remote
`structural-fallback-source-v1`, manifest SHA256
`90c1faee28214c5a41a9071d93b0ece6242028728c2b550653da167c123540f7`;
new prospective scaling protocol SHA256
`0f3cd4a7f3280932f3faacc31ccf26d3dbba43ed7b133bf8eacc8ee8722a3d74`.
Six tests pass using the actual failed-root and valid-root fixtures, checking
exact fallback behavior, sample/order invariance, unchanged valid inference,
propagation of other errors, and incomplete/duplicate pair rejection. Parity
job 10858908 completes in eleven seconds: all 623 original calibration labels
and all 585 original development labels remain unchanged, without fallback.
The new rule emits all 1,552 small-cohort pairs and scores 1,510 correct
(97.2938%), with fallback only on cohort 000. These small cohorts are now
explicitly development diagnostics for the changed policy; they are not
independent qualification. The original calibration scalar remains 0.7744
because all calibration inference is unchanged. It still has zero individual
resolution and does not pass an ECE gate merely because labels improve.

The corresponding inference handoff v4 (unlaunched v3 packaging draft retained)
adds the narrow policy and explicitly verifies the full fallback source
closure. Its protocol SHA256 is
`b71a7d3dc1970318b898c0e1187dd691517d7d7ca696412b0ff6e9f016d886ad`;
manifest SHA256 is
`07b62a05619b3d941ea8d0cbb46f8cab2b1146fd3198aac12e91a04dfc4b7613`.
The 263-person full-roster contract, map crop, input phase/missingness policy,
model, features and CSV rules are unchanged. The two original 229-person
cohorts remain genetically unrun at 20:53 UTC; source and launch-contract
review plus independent new-policy diagnostic rescoring are required before
launch. Each retains eight CPUs, 90 minutes, 44 GiB and the 21:50 UTC actual
start cutoff. No final data has been accessed by this agent.


At 20:55 UTC the parent completed independent source, launch, parity and raw
row review and authorized the original two fresh 229-person cohorts under the
new policy. Array 10859039 (indices 8 and 9, concurrency two) uses the unchanged
8 CPU / 90 minute / 44 GiB per-cohort cap. The immutable launch record SHA256 is
`68a2a5de391e52cf4c75ef218ad51e327d4a969d373c67b3683a82e2125bc3ba`.
It binds the parent's review and original strict failure. The protocol's inherited
`model_or_method_changes:false` field is explicitly interpreted by both review
and launch as `method_change:true, fitted_model_change:false`: adding an
exception-only fallback changes the inference policy, not the fitted forest.

Fixed assessment job 10859075 waits for both array jobs and requires both full
rosters, 52,212 predictions per method, raw output hashes and Slurm accounting.
Its new source manifest SHA256 is
`7946ccb15c64fb72a4257bdc9235bc94d02c2e1e549b6b90e0b8b54a1c12712a`.
The pragmatic comparison is reported for each cohort and pooled: exact accuracy
at least 0.70, macro recall at least 0.60, weighted Jaccard above the matched
baseline, complete coverage and all resource caps. Original strict numerical
and ECE diagnostics remain separate and cannot reverse the original strict-v3
failure. Export job 10859110 follows assessment and retains all raw truth,
rosters, structural/pair/baseline predictions, tree provenance, logs, manifests
and missing-file records.

At 21:07 UTC the parent conveyed the user's explicit deadline-first instruction
to submit remaining entries regardless of scientific gate status. This changes
the submission decision, not the thresholds or evidence: the same frozen policy
and fresh assessments continue, with failures disclosed and no additional
synthetic tuning. The parent owns final 263-person inference and submission;
this agent has not accessed final input data or predictions.


At 21:10 UTC the parent authorized this agent to read frozen final output
metadata and predictions solely for official packaging. No final VCF access,
new inference or tuning was authorized or performed. Parent-owned inference
job 10859028 completed and supplied all 263 people / 34,453 pairs. Packaging
script `package_task10_final_v1.py` verifies the original preflight roster and
input identity, frozen source manifest and validator, exact raw CSV equality
to the frozen prediction rows, and before/after input hashes. It appends only
`meioses_count`, using frozen HUMAN_LABELS `up + down`, independently equal to
`relation_meioses(10)`. It preserves every original ID, relation and row order.

The separate official four-column output is
`final_runs/task10-root-20261120-2028/task10_final_candidate.csv`, SHA256
`ee9bace7bfc086281be0b243c5b440db9542277466b4722a14cb32523ef059be`.
Frozen validation confirms all 34,453 exact full-roster unordered pairs, 263
unique individuals and matching relationship/meiosis fields. Packaging
provenance SHA256 is
`1643526d98cb9d270f53f493c66271500cff94164a718049123f99a6e16bb931`;
its local mirror is `task10_final_packaging_provenance.json`. The script passes
Ruff. Original three-column candidate and prediction files remain unchanged.
This agent performs no upload or submission, and packaging does not change
the original strict qualification failure. The parent independently reviews
and owns the one-shot submission decision.


The two fresh 229-person cohorts and all assessment/export jobs completed.
Each cohort contains all 26,106 pairs, and neither used fallback. Structural
accuracy is 17,598/26,106 (67.4098%) and 17,786/26,106 (68.1299%); pooled
35,384/52,212 (67.7699%), macro recall 77.3713%, minimum-label recall 15.0794%,
and weighted Jaccard 0.5244112. Matched pair-only accuracy is 40.4543% and
weighted Jaccard 0.2890453; original baseline accuracy is 32.9369% and weighted
Jaccard 0.2320804. All complete-pair, macro, relative-WJ and resource gates
pass; the prospective practical exact-accuracy >=0.70 gate fails in both
cohorts and pooled. Strict numerical failures are also retained. No tuning or
resubmission follows these outcomes.

Array tasks complete in 28m56s and 28m42s, with peak RSS 7,005,252 and
7,145,192 KiB, below the declared caps. The fixed assessment completes in
31 seconds and raw export in 11 seconds. Export manifest SHA256 is
`19b57b168b8f05165e815ef4ea56544a2177d4ff1e3ad0b70cb6d6940b0ebdeb`: all
29 copied files match, no files are missing. All raw truth, roster, three
methods' predictions, score inputs, provenance, source manifests and accounting
are mirrored under `structural_fallback_scaling_evidence/`. The parent receives
these for independent audit.

The parent reports Task10 submitted at 21:24 UTC as `syn78780148v1`, with UI
success, under the user's deadline-first instruction. This agent performed no
upload or UI action. Submission does not imply strict or practical qualification
passed. A compact record is `task10_final_evidence_summary.json`.


Correction to the earlier UI submission-success record: the parent subsequently
received server rejection for submission 9785356 in queue 9619674 (email
`1a1fad9451a8e239`, 21:25 UTC). The entry was invalid and not scored: the actual
queue accepts 12 human relationship labels, while the generic local validator
allowed 14. This is a schema-validation failure; earlier statements that the
CSV passed local validation do not establish server acceptance. The invalid
CSV is preserved unchanged.

At the parent's explicit instruction, `correct_task10_server_schema_v1.py`
implements only the server allowlist correction. It maps unsupported
`full_1st_cousin_thrice_removed` to `full_2nd_cousin_once_removed` (both seven
meioses), and `avuncular_twice_removed` to `full_1st_cousin_once_removed` (both
five). Each is the unique allowed same-meiosis label. These are different
pedigree relationships; the adapter merges labels for the queue schema and
does not claim biological equivalence. Exactly 3,040 rows use the first mapping;
zero rows use the second. All 34,453 pairs, 263 IDs, row order, meiosis counts
and other labels are unchanged. A separate hard-coded 12-label/count checker
passes independently of the adapter. No model, threshold, final-truth or
scientific score feedback is used.

The corrected candidate SHA256 is
`42ac7e0de7953d7f79f7cd6ad3fcc356c963f43fae6ee636646797fb90639eb8`.
It and correction provenance/independent check are copied byte-identically
to the parent's `/private/tmp/ghist-final-submissions-20261120/task10/` and
`experiments/final_submission_2026_11_20/task10/`. The parent owns eligibility
checks, UI and any retry. This agent has not resubmitted. Original scientific
and validation failures remain unchanged; compact correction record is
`task10_final_schema_correction_summary.json`.
