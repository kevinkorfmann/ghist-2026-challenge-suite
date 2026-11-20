# Task 11 independent varied-pedigree development

This is ongoing development, with no final challenge data or submissions.
The frozen v3 run failed its scientific qualification criteria and is not ready
for a final challenge run or submission.
The previous sibling-only examples are not reused for fitting or confirmation.

The official Task 11 description specifies eight unrelated dog families,
possible within-family inbreeding, and one population of originating individuals.
It does not identify an empirical breed founder panel. This study therefore uses
independent coalescent founder haplotypes with LD and the existing dog species map,
with preregistered population-size and missingness stresses. This is an explicit
model of the unspecified population assumptions, not evidence of a breed match.

The new candidate preserves pair-model predictions. Bonsai is not allowed to
replace supported labels. Complete-link, strongly supported full-sibling groups
can optionally pool relationship support to other sibling groups, but only for
weak original calls and only if the separate calibration cohorts show improved
balanced accuracy without introducing errors among previously correct calls.
All numerical rules are fixed in the prospective protocol. The policy is frozen
before any new confirmation cohort can run.

The 16 training, four calibration and eight confirmation cohorts have disjoint
seeds, founder draws and jointly simulated pedigrees. Each cohort contains all
eight accepted relationship labels, eight families with variable branch sizes,
and unrelated, first-cousin and full-sibling ancestral mating. The pipeline
simulates phased genotypes, introduces fixed MCAR missingness, drops incomplete
sites, executes real hap-IBD, corrected diploid conversion, HBD masking, family
clustering and candidate classification. Every expected true pair remains in the
scoring denominator; extra cross-family output is separately counted.

The prospective confirmation gates are overall accuracy at least 0.80, balanced
accuracy at least 0.75, each-label recall at least 0.60, and a cohort-bootstrap
lower 95% balanced-accuracy bound at least 0.70. At least 95% of whole families
must be recovered exactly, with no missing samples or true pairs and no extra
cross-family output pairs. Fixed-bin expected calibration error must be at most
0.10. These are declared minimum utility targets, not official challenge rules.
They will not be lowered after observing confirmation outcomes.

Implementation: `scripts/develop_task11_pedigree.py`. Evidence and the copied
prospective protocol are retained in `experiments/task11_development_2026_11_20/`.
All substantive simulations run in the isolated Betty experiment directory
`/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120`.
No shared production source, existing evidence, or frozen release is modified.

The v3 protocol and source snapshot were frozen before confirmation. The
training-only clustering comparison selected average linkage at the fixed 10 cM
affinity threshold: it recovered 128/128 training families, compared with
125/128 for the previous spectral method. The candidate was then refit on 2,848
training pairs. Four separate calibration cohorts gave 680/735 correct calls,
85.36% balanced accuracy and ECE 0.02768. Temperature 0.75 was selected;
sibling pooling did not meet its calibration condition and remained disabled.

All eight prospective confirmation cohorts completed the genotype-to-detector
pipeline. The retained evidence is in `evidence_v3/` inside the experiment
directory. There were 1,443 true pairs, 1,223 correct predictions (84.75%),
73.82% balanced accuracy, and a cohort-bootstrap 95% interval of 68.40–79.15%.
Avuncular-once-removed recall was 52/116 (44.83%). Balanced accuracy, minimum
label recall and the bootstrap lower bound failed their frozen gates. All
64 families were recovered exactly, every required pair was present, no
cross-family extra pair was emitted, and ECE was 0.03801. These successful
coverage checks do not override the failed scientific gates.

After failure, descriptive diagnostics found once-removed recall of 23/42 in
the four simple cohorts and 29/74 in the four complex cohorts. Recall was
5/34 with unrelated ancestral matings, 18/46 with first-cousin background and
29/36 with sibling-mating background. Those strata also corresponded to
maximum detected endpoint HBD below 5%, 5–15% and at least 15%, respectively;
they are overlapping descriptions, not independent effects. The major
confusions were half-avuncular (28), full first cousin (17), and half first
cousin (11). These cohorts are now development evidence and cannot be reused
as fresh confirmation for a revised candidate.

The separate `scripts/prepare_task11_candidate.py` adapter binds the candidate,
protocol, source snapshot and model hashes, preserves the VCF genotype and
phase fields while assigning reversible internal sample aliases, invokes the
frozen detector/features/clustering, restores original sample IDs, and exports
`task11_candidate.csv`. It validates all within-family pairs with the existing
frozen Task 11 validator. Missing or failed qualification evidence requires
explicit development-only mode. Ten adapter boundary tests pass, and Ruff
passes. One synthetic full-pipeline rerun reproduced all 183 labels with a
maximum probability difference of 1.67e-16. Separately, 160 ID/order
permutations preserved the eight confirmation cluster partitions; these are
cluster-only metamorphic checks, not 160 classification reruns or independent
biological cohorts. The adapter does not upload or authorize any final data.

The failed candidate remains immutable. A proposed next development step is
to measure how IBD is distributed over the two parental haplotypes, masking
IBD2 and HBD. This follows the rationale of the chromosome-wise haplotype
concentration statistic in [PONDEROSA, equation 2](https://doi.org/10.1093/genetics/iyaf094).
Any resulting improvement will require new, prospectively recorded
confirmation seeds and the same numerical gates.

The minimum phase extension is now frozen as source snapshot v4 (method v3),
with prospective protocol SHA256
`da6f44548fdf8e02844221b780bc8585d7b2cddc1ba369fa2000de364acff4e5`.
At each fixed HBD cutoff of 2, 5 and 10 cM, it adds the sorted two endpoint
concentrations, retained cM denominators, and retained/pre-HBD exposure
fractions. True IBD2 uses two distinct homolog matches and is always masked;
overlapping copies of the same edge cannot manufacture IBD2 or extra mass.
Zero exposure is encoded explicitly as zero concentration, denominator and
retained fraction. Concentration is computed per chromosome before summing,
so chromosome-wise homolog orientation is irrelevant. The 37 development
tests include endpoint swaps, independent chromosome homolog flips, duplicate
and overlap invariance, HBD masking, genuine IBD2, one-homolog fans, zero
exposure and map-bound failures. No detector, clustering or classifier
hyperparameter changed.

The 16 original training cohorts and four calibration cohorts completed this
frozen replay. Calibration decreased to 666/735 correct calls and 81.44%
balanced accuracy; once-removed recall decreased from 32/39 to 25/39.
On the eight exposed development cohorts, the phase model improved to
1,254/1,443 correct calls and 76.85% balanced accuracy, but once-removed recall
was only 58/116 (50%) and half-avuncular recall 42/78 (53.85%). Therefore no
fresh confirmation was run. The matched ablation verified equality of all
original features, and the eight exposed cohorts were not used for fitting.
The separate adapter v2 accepts the explicitly bound old and phase protocol
versions; its 11 boundary tests and Ruff pass.

A single further phase hypothesis is frozen as source snapshot v5 (method
v4), with protocol SHA256
`ca55f4dd692e0c1aa41704dcb3901b30b129437b1e0905211ceb84f00bfd0d1c`.
It requires at least 10 cM of merged support on an individual homolog edge
before that edge contributes phase exposure. The IBD2 mask still uses all
detected segments, including shorter opposite-homolog edges. This is motivated
by a separate Task 10 comparison on identical genomes: short detected sharing
altered haplotype fractions substantially, whereas total sharing above 10 cM
was much more stable. It is a hypothesis requiring dog-specific validation,
not evidence of transfer. All 18 phase features, the three HBD cutoffs, explicit
zero-exposure behavior, detector, base features, clustering, classifier and
scientific gates otherwise remain fixed. All 38 development tests and Ruff
pass, including short-edge exclusion and masking by a short genuine IBD2 edge.

The v5 development replay reuses immutable v4 detector outputs and preserves
all 131 base features byte-for-byte, recomputing only the 18 phase features.
Its frozen helper verifies the predecessor's inputs, source hashes, all bound
output hashes and the unchanged detector/model function ASTs. Each new result
records this lineage and explicitly says that the detector did not execute
again. The same separate training and calibration cohorts are used; the eight
exposed cohorts are reserved for a matched development check. Fresh
confirmation seeds 2026120400–2026120407 remain unused until the unchanged
development criteria pass. Any new qualification outputs will use their own
`qualification_v5/` directory to preserve earlier evidence.

The v5 matched development check completed with 1,254/1,443 correct calls,
76.76% balanced accuracy and ECE 0.02125. Once-removed recall increased to
66/116 (56.90%), but half-avuncular recall was 43/78 (55.13%). Both remain
below the fixed 60% minimum; no fresh confirmation was run. Calibration was
673/735 correct with 83.82% balanced accuracy. The adapter's explicit schema
support now covers this frozen variant; all 12 boundary tests and Ruff pass.
No new full adapter smoke or scientific qualification is claimed.

A separate fixed information diagnostic now measures union sharing across
predicted sibling groups and IBD^(011): regions where two predicted siblings
both share with a relative but have no detected IBD with each other. This is
motivated by [DRUID](https://doi.org/10.1016/j.ajhg.2018.05.008), whose human
relationship-degree results do not validate these eight dog labels in inbred
families. Training sibling groups come from four cohort-level out-of-fold
models. An entire ambiguous non-clique component falls back to singletons.
No truth grouping, sample ID or pedigree property enters its 36 features.
Merged homolog edges of at least 10 cM contribute sharing; all detected IBD
between the putative siblings excludes IBD0, and HBD from all three individuals
is masked at each existing cutoff. Zero exposure is explicit. One diagnostic
model retains the same ExtraTrees settings and calibration grid, and all
target pairs are scored against the matched frozen v5 baseline. The exposed
eight cohorts remain development-only. The prospective plan and implementation
are retained as `context_diagnostic_plan_v1.json` and
`diagnose_family_context.py`; this diagnostic is not a replacement candidate.

The context diagnostic completed with all eight interval/invariance tests and
Ruff passing. On the exposed cohorts it increased accuracy to 1,312/1,443
(90.92%) and balanced accuracy to 82.91%. Once-removed recall reached 71/116
(61.21%), but half-avuncular recall decreased to 39/78 (50%). Calibration
remained 673/735 correct, with once-removed recall 20/39 and half-avuncular
recall 31/56. It corrected 97 exposed errors and introduced 39, including two
previously correct predictions with at least 90% base-model support. These
results show useful observed family information but do not justify promoting
the model. All raw scorer inputs, grouping audits, features and the fixed plan
are retained in `context_diagnostic_v1_compact.tar.gz`; model weights remain
hash-bound in the isolated remote directory. No confirmation was run.

The next bounded experiment expands independent training from 16 to 64 whole
cohorts, adding seeds 2026121000–2026121047. Each accepted label has 64 focal
families rather than relying on duplicated correlated pair rows. The combined
training set has 32 simple and 32 mixed-generation cohorts. The fixed Ne and
missingness grid, 2–4 member branches and three ancestral-mating backgrounds
remain unchanged; each new focal label receives 16 pedigrees in each mating
background. Four cohort folds each contain eight simple and eight mixed
cohorts. The original four calibration and eight exposed cohorts remain
development data, and all reserved confirmation seeds remain unused.

`expanded_training_plan_v2.json` and `expanded_training_bundle_v2.tar.gz`
freeze the generation protocol, source closure, model settings, split rules
and resource limits before generation. All 48 new cohorts must complete and
match their registered identities; there are no outcome-based exclusions.
The fixed pair model and context model will be refit on the identical 64
cohorts, with the same calibration temperature grid. The generation array is
capped at four concurrent tasks, each using four CPUs and 5632 MB per CPU on
Betty; every new cohort executes genotype generation, hap-IBD detection,
family clustering and feature extraction. This is still development, not a
qualification claim.

An additional eight-cohort sparse-family stress design is frozen before the
expanded fit finishes. Seeds 2026122000–2026122007 are disjoint from all prior,
expanded-training and reserved-confirmation streams. Four simple and four
mixed-generation cohorts use focal branch sizes 1–4, with all eight focal
labels and the original population, missingness and ancestral-mating settings.
`sparse_context_stress_design_v1.json` and its explicit launch-guard addendum
require the expanded context model to pass all nine original gates on the
exposed development cohorts, plus calibration minimum label recall 60% and
ECE at most 0.10. There is no calibration bootstrap requirement. The stress
run will bind the completed expanded fit's exact models and temperatures;
it permits no refits or outcome exclusions. Stress must pass all original
gates before the reserved confirmation seeds can be used. Neither stress nor
fresh confirmation has been launched.

The sparse design's pedigree metadata passes 513 independent ancestry-label,
pair-coverage and seed checks, with Ruff passing. Its registered seeds imply
1,204 pairs across 64 families; all eight labels are represented, including
98 once-removed and 62 half-avuncular pairs. This preview uses no genotype,
detector or classification outcomes and does not count as validation.

All 48 added full-pipeline cohorts completed, with all 384 new families
recovered exactly. The fixed expanded model fit is in progress. The separate
v5 adapter smoke also completed: it reran the synthetic VCF through detection,
clustering, classification and validated CSV export with 183 pairs, zero label
changes and maximum probability difference 1.11e-16. This reuses one exposed
biological cohort and verifies the adapter, not independent scientific
qualification. The result is `adapter_smoke_report_v3_development.json`.

The expanded 64-cohort fit failed. Its context model scored 1,231/1,443 on
exposed development, with balanced accuracy 0.71860, half-avuncular recall
19/78 and bootstrap lower bound 0.62648. Calibration minimum recall was 0.50.
The matched pair model also regressed, to balanced accuracy 0.65128. The raw
launch audit denied stress execution; neither the sparse design nor any fresh
confirmation was run. Full raw metadata, predictions, context features and
hashes are preserved in `expanded_context_development_v2_compact.tar.gz`.
The context inference wrapper separately reproduces 203 retained calibration
predictions per model with zero label changes and differences below 2.3e-16;
this parity check does not change the failed scientific result.

An audit of all 76 retained cohorts and 608 families verified a design flaw:
within each split and focal label, population size uniquely determined the
ancestral mating kind. All 24 label-by-population cells rotated between the
original and added training blocks because the new local population index
differed from the pedigree index. The mapping is verified in
`nuisance_schedule_audit_v1.json`; its causal contribution to the performance
regression remains a hypothesis. Earlier independent genomes did not remove
this shared nuisance correlation.

The corrected prospective design freezes 108 new training, 54 calibration and
54 development cohorts. Each block crosses three population sizes, three
missingness rates, three ancestry shifts and two topology types. Training has
two independent blocks; each label therefore receives every ancestral mating
kind within every population, missingness and topology cell. Branch sizes are
1–4 throughout. Explicit ancestry and third-branch assignments replace the
cohort-index coupling; optional training third branches are balanced between
replicates. Four held-out training folds each contain 27 cohorts, every
population/missingness/ancestry-shift cell once, and 13/14 simple/mixed
topologies. All prior data are excluded from fitting and calibration.

`crossed_development_plan_v1.json` and `crossed_development_bundle_v1.tar.gz`
bind the schedule, new source, unchanged v5 detector/features and unchanged
context classifier. Nine design/mechanics tests and Ruff pass. All 216 outcomes
must be retained; thresholds are unchanged. New stress seeds 2026126000–6007
and confirmation seeds 2026127000–7007 are reserved before outcomes, with
their full configurations recorded. Older reserved seeds are not the corrected
confirmation set. Generation is capped at four four-CPU tasks; hard start
limits are 21:40 UTC for generation and 21:50 UTC for fitting, with a 22:15
development decision and 22:45 handoff target. Root independently verified
all source bindings, schedule crossings, folds and master-seed separation at
19:42:37 UTC. Generation job 10856905, dependent fixed fit 10856907 and raw
audit 10856920 then launched. These remain synthetic finite-template
assumptions, not a claim of universal pedigree coverage.

The separate `crossed_future_operational_bundle_v1.tar.gz` freezes a guarded
stress/confirmation runner and truth-free context CSV adapter before any
reserved outcome. Twenty-one fail-closed gate/binding tests and Ruff pass.
Source-only initialization verifies the complete v5 snapshot before scientific
imports; candidate models will be bound only if the unchanged development
launch gates pass. All eight configured stress cohorts must pass before the
eight reserved confirmation cohorts can start. The adapter preserves VCF
phase/missingness handling, exports all inferred within-family pairs and uses
the existing CSV validator. No candidate, reserved run, final-data access or
submission has been created by these preparations.

At 20:21 UTC the user accepted reasonable imperfect candidates for submission.
This changes the submission decision, not the scientific gates. A separate
pragmatic package now binds the existing 16-training-cohort context model,
without refitting or claiming qualification. On the eight exposed cohorts it
has 1,312/1,443 correct (90.92%), balanced accuracy 82.91%, weighted Jaccard
0.85669 and ECE 0.02324. Its matched pair baseline has 1,254/1,443 correct,
balanced accuracy 76.76% and weighted Jaccard 0.80105. Always predicting
siblings on those same pairs gives 50.52% accuracy and weighted Jaccard
0.25522. These are synthetic development comparisons, not expected final
challenge scores. Half-avuncular recall remains 39/78, and calibration also
fails the half-avuncular and once-removed minimum-recall gates. The old
generator's confounding and sibling-rich structure remain material limitations.

`pragmatic_context_candidate_v1.json` has SHA256
`6ccfda5f8133b5dee5265f35ac3a328fc4c35cfa88af076e11b9873af781ea15`.
Its distinct schema explicitly records failed strict qualification and cannot
pass the frozen future-cohort runner. Fifteen identity/hash/scope tests and
Ruff pass. The truth-free CSV adapter was rerun through the complete detector
on one retained calibration cohort after simultaneously changing all sample
aliases and column order. Job 10858308 completed in 68 seconds: all 203 pairs,
all eight family partitions and both models' labels were unchanged; maximum
probability differences were 2.22e-16 and 1.67e-16. Historical context
predictions were also reproduced. The existing CSV validator passed. This is
one operational classification rerun, not new biological validation.

The compact handoff bundle `pragmatic_context_handoff_bundle_v1.tar.gz` has
SHA256 `3230983b1105b181e0be74d6e91cb008d7b736eebc8d4408caf091dad2f77e64`.
Root received the fixed adapter, complete scores, limitations and a four-CPU,
22-GiB, 45-minute inference template. Root owns final data and the one allowed
final submission; this agent has neither accessed final data nor submitted.
The candidate is an operational fallback, not a commitment to final selection.

Before corrected outcomes, a fixed comparison was frozen to evaluate the old
base/context models on every one of the 54 corrected development cohorts.
It reuses observed detector outputs, fits nothing, retains every outcome and
uses no reserved seeds. The plan SHA256 is
`de370d0e5aeb5246ec7d34ec19864a1ba7a8505669ceb5ca2ff325f80d6cde7d`;
job 10858377 depends on completion of all 216 generation tasks and is bounded
to four CPUs and 15 minutes. This supplies a matched distribution comparison
for root's single final choice. Selection on these results makes them
development evidence; it does not constitute independent confirmation.

Before that comparison executed, self-review found that its auxiliary Jaccard
calculation needed to count extra inferred-family pairs as false positives.
Pending job 10858377 was cancelled without running; its files are preserved.
Replacement job 10858576 uses the same 54 cohorts and fixed models, with raw
truth/prediction scoring and four targeted missing/extra/duplicate tests. The
v2 comparison plan SHA256 is
`ed744cc5c2cf3a0198b8534748d216784519d3a1f70e46436ebec958ba9acbcc`.
All strict jobs, criteria and reserved cohorts remain unchanged.

Root's first final-data preparation failed before classification because a
detected interval exceeded the frozen map's physical endpoint. The failed run
is preserved. A separate external staging step now crops VCF records to the
exact inclusive physical domain of that same map. It preserves every retained
header, genotype, phase and FORMAT byte, rejects unmapped chromosomes, checks
at least 50 retained raw records per map chromosome and records all exclusions
by chromosome. It performs no cM extrapolation, imputation, model change or
scientific-source change. Subsequent missingness/phase preparation remains the
unchanged frozen adapter's responsibility.

Seven staging tests and Ruff pass. Synthetic job 10858770 completed in 62
seconds. Cropping retained calibration_00 excluded nothing; adding one copied
out-of-map genotype row per chromosome and cropping again removed exactly 38
records and yielded the same deterministic gzip bytes as the identity crop.
One full unchanged-adapter rerun reproduced all 203 base/context labels and
all eight families, with maximum probability differences below 2.3e-16; the
CSV validator passed. This is an operational domain-boundary check, not new
biological qualification. Root received `pragmatic_context_map_crop_v1.sbatch`
and `map_bounds_handoff_bundle_v1.tar.gz`, bundle SHA256
`abbc43be2d1a13131442e7c451fce00dd85d1ef89acdf522418a62027952d8e5`.
This agent did not access or modify final genotype data or perform a retry.

### Deadline-first fixed comparison, 21:12 UTC

The user authorized submitting all remaining tasks even if strict validation is incomplete. Root retains sole control of final data, single-entry selection, and upload. This authorization does not change the recorded scientific gates. Before any corrected classifier outcome was available, root fixed selection on the same 54 corrected synthetic development cohorts across old base/context and corrected base/context: complete valid outcomes, accuracy at least .70, macro recall at least .60, minimum class recall at least .30, family recovery at least .95, then highest official truth-weighted Jaccard (ties within 1e-12 prefer old context, new context, old base, new base). Under the later deadline-first instruction, root may choose the highest valid Jaccard if none passes those practical thresholds. No final-score feedback, per-pair mixing, refitting, or threshold tuning is allowed. The timing record confirms the corrected fit and comparison had not produced any classifier outputs when the practical rule was recorded.

All 216 crossed cohorts completed successfully, with 1728/1728 family partitions exactly recovered and explicit ancestry assignments checked. The old-model fixed comparison completed on every one of the 54 development cohorts and all 7305 true within-family pairs. Old base scored 5690/7305 (.778919), macro recall .631035, minimum recall .315789, official weighted Jaccard .695161. Old context scored 6004/7305 (.821903), macro recall .693807, minimum recall .332632 (half-avuncular), official weighted Jaccard .748240. Neither has missing or extra pairs. These broader, sparser, independently crossed results supersede historical development scores as the basis for the fixed final comparison; they are development selection evidence, not independent post-selection confirmation.

The full old-model comparison, raw truth, rosters, partitions, and bindings are preserved in `pragmatic_crossed_comparison_v2_compact.tar.gz` (SHA256 `9b021cba36d138fde4217084d9e62ade8ede2ed5847b5c937a28bad0065af667`). Corrected fitting and raw auditing continue unchanged. A dependent Slurm job (10859347) will bind both corrected model roles without selecting either for submission and will run one complete retained-synthetic alias/order parity check for both models. Reserved stress and confirmation cohorts remain unused.

### Corrected fixed-model handoff, 21:21 UTC

The corrected context model is recommended for root's one-time selection under the prospective practical policy. All four candidates have complete valid predictions on the identical 54 development cohorts (7305 pairs and 432 exactly recovered families). The fixed comparison is:

| Model | Correct pairs | Accuracy | Macro recall | Minimum recall | Official weighted Jaccard |
|---|---:|---:|---:|---:|---:|
| Old base | 5690/7305 | .778919 | .631035 | .315789 | .695161044 |
| Old context | 6004/7305 | .821903 | .693807 | .332632 | .748240362 |
| Corrected base | 6131/7305 | .839288 | .718837 | .421053 | .760241094 |
| Corrected context | 6535/7305 | .894593 | .808057 | .570526 | .835006231 |

Corrected context recall is .999718 for siblings, .934066 for half-siblings, .975325 for avuncular, .570526 for half-avuncular, .743777 for avuncular once removed, .794411 for avuncular twice removed, .763214 for full first cousins, and .683417 for half first cousins. Its ECE is .007592 and cohort-bootstrap macro-recall 95% interval is [.782564, .833296]. It corrects 565 errors from the matched corrected base while introducing 161 errors, including one change from a previously correct prediction with base support at least .90.

Strict qualification remains false. The half-avuncular development recall fails the original .60 requirement; its calibration recall is 170/405 (.419753), also below the launch guard. No reserved stress or confirmation was launched, and root explicitly stopped additional development/confirmation for the deadline-first handoff. The 54-cohort comparison is development selection evidence. Sparse branch sampling, fully crossed nuisance settings, LD, missingness, and inbreeding are represented, but finite templates, true phase, MCAR missingness, and coalescent founder assumptions remain transfer limitations.

The full corrected raw archive is `crossed_development_compact_v1.tar.gz` (SHA256 `f282a52aec52f29bc99262f344766e3b394d7f5d60db982b92ada9265b145339`). The context candidate is `pragmatic_corrected_context_candidate_v1.json` (SHA256 `f2c1391a5866567a3bc10423bea25ecd537b86f9d3cf819fcd72798bb2abd833`), binding 679 inputs including exact source, models, calibration choice, raw reports, and retained failed scientific audit. The context model SHA256 is `ec8becd738f28ec3b494a408887ee1c16303d69ac872c39ceae46f3488f77ef5`, with temperature .75. The base model is bound separately with temperature 1.0. Packaging does not select the final submission; root records that decision independently.

The unchanged corrected CSV adapter passed one complete synthetic genotype-to-detector-to-clusters-to-model-to-CSV rerun after simultaneously randomizing sample aliases and genotype column order. Both models preserved all 129 pair labels, with maximum probability differences 2.22e-16 (base) and 1.11e-16 (context); all eight family partitions were unchanged and the CSV validator passed for all 46 samples. This is operational parity, not additional biological validation. Slurm job 10859347 completed in 90 seconds, and the compact export 10859432 completed in five seconds.

The 26-file operational bundle `pragmatic_corrected_operational_bundle_v1.tar.gz` has SHA256 `6b41fabb303c43bd54972a6cd06663031fab9a8c3b389edcef3bce5c688408e2`. It includes both candidate manifests, the exact source/initializer, raw synthetic predictions, official three-column CSV, and audit/validation evidence. The root-only wrapper `pragmatic_corrected_inference_v1.sbatch` (SHA256 `3f37a40ec61b6725a8e215655d70ceba9b9f27fb5dafac423d021f6f51253daf`) requests four CPUs, 5632 MB per CPU, and 45 minutes. Its positional arguments are an already map-bounded VCF, a new output directory, a fixed role (`context` recommended), and the exact candidate SHA. Root can reuse the independently verified map-bounded input from its earlier preparation. This agent did not access final genotypes or predictions, select or upload a final entry, or alter the original failed final run.
