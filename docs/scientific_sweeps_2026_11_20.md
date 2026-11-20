# Scientific validation of Challenges 5–8 — November 20, 2026

No public challenge answers were accessed for this work. Final-data inference and
submission remain deferred. This review distinguishes reproducible software,
retrospective simulation reanalysis, and independent biological evaluation. The
archived freeze files are historical records and are not scientific certification.

## Material findings and corrections

The old simulation recovery pipeline chose the fine-localization region using
the simulated truth ±750 kb. Both the largest-gap and split-flank linkage
estimators therefore had information unavailable in a real prediction. This
contaminated hard-sweep ensemble validation and soft-sweep method selection.
The revised `ghist_single_sweep.recovery.localize_from_scan` chooses the region
from the full composite scan and has no truth-position argument. Truth is used
only after prediction to score errors. Records identify the search rule and
state that inference did not use truth. New calibration and freeze entry points
reject old records without this provenance and reject duplicate simulation seeds.
The fixed composite-region rule was chosen before inspecting revised errors; it
is an explicit reproducible candidate, not validation of the manually selected
historical testing regions.

The independent RAiSD comparison also paired lexically sorted simulation files
with numerically indexed reference files. For example, `task_10` could receive
`task_2`'s reference estimate. It now joins by filename and checks the generating
truth identity. Simulated VCF export now writes 1-based positions explicitly and
converts the RAiSD result back to 0-based positions. The simulation CLI records
its mutation rate, demographic model, and stdpopsim version; new runs default to
exact SLiM scaling (Q=1). Its default rate is
now the documented testing rate, 1.44e-8, and can be explicitly overridden.

A requested recombination map could silently become a flat map when downloading
failed or concurrent downloads temporarily left the chromosome file missing.
The scientific runner and sweep simulator now reject this fallback. The first
new BGS smoke attempt exposed this failure and a missing RAiSD GSL runtime path;
its results were discarded. The rerun used the actual named map and the required
GSL library. Numerical population sizes and times are passed directly to the
engine; no compensating pre-scaling is performed. stdpopsim itself warns that
scaling is approximate, especially with selection. See the
[stdpopsim SLiM implementation](https://popsim-consortium.github.io/stdpopsim-docs/stable/_modules/stdpopsim/slim_engine.html)
and [selective-sweep API](https://popsim-consortium.github.io/stdpopsim-docs/stable/api.html).

Challenge 7's exporter still transferred Challenge 5 recovery, despite BGS and a
different ten-estimate method. It now requires `--matched-recovery` from the
matching bottleneck/growth plus BGS scenario and the exact ten-estimate method;
the archived transfer cannot pass. Matching declarations are necessary metadata,
not proof of adequate simulation coverage. Challenge 8's existing independent
null gate remains necessary; its statistical exchangeability assumption and power
must be tested on the intended chromosome and nuisance-parameter range.

## Hard- and soft-sweep reanalysis

All 192 archived exact-scale hard-sweep simulations and all 48 standing-variation
soft-sweep simulations were recomputed with truth-blind localization. Whole
simulation seeds remain the split units. These are previously evaluated
simulations, so this is a retrospective bug-fix reanalysis, not a fresh untouched
test. The stored `untouched_evaluation` key remains for schema compatibility;
`evaluation_status` explicitly records this limitation. No revised method was
frozen from these results.

| Challenge | Revised evidence | Scientific interpretation |
|---|---|---|
| 5 — Single hard sweep | 56/192 ensemble estimates changed. The final 64 seed-heldout simulations have median error 25.38 kb, 90th-percentile error 118.56 kb, and maximum 368.74 kb. Of these, 11 pass the pre-existing 50 kb estimator-consensus rule, and all 11 lie inside the 60 kb calibrated half-width. | This supports limited conditional localization in the tested 5 Mb simulation domain. Eleven conditional successes do not establish general coverage, and the region is much smaller than a full chromosome. |
| 6 — Soft sweep | For the historically selected gap estimator, the final 24 simulations still have 13/24 coverage inside 35 kb, but the 90th-percentile error increases from about 678 kb to 2,023 kb after removing the oracle region. Repeating the original calibration procedure selects composite instead: median heldout error 39.49 kb, 90th-percentile error 2,187 kb, and 16/24 coverage at its selected 75 kb half-width. | The broad scientific gate fails. Narrow signal-conditioned summaries cannot replace this failure across weak and incomplete sweeps. The revised calibration choice is diagnostic and is not a promoted prediction method. |

Full records and generator seeds are in
[Task 5 reanalysis](../experiments/scientific_2026_11_20/sweeps/task5_blind_recovery.json)
and [Task 6 reanalysis](../experiments/scientific_2026_11_20/sweeps/task6_blind_recovery.json).

## Independent biological background-selection benchmark

A separate runner generates new SLiM simulations on chromosome 7: 36–37 Mb,
HapMapII_GRCh38 recombination, actual HAVANA104 exon intervals, deleterious
Gamma_K17 mutations, 25 diploid samples, and mutation rate 1.44e-8. The arbitrary
history has ancestral size 6,000, size 1,200 beginning 1,200 generations ago,
and recent growth over 300 generations to 12,000. Positive controls have either
one or two exonic selected sites, selection coefficient 0.04, origin 650/750
generations ago, and conditioned end frequency at least 0.95. These are scenario
matches, not estimates of the unknown challenge-generating history. Samples,
seeds, exon hash, model settings, and report hashes are recorded.

Before examining positive evaluation results, the detector was fixed at alpha
0.05, separation 100 kb, half-width 50 kb, and maximum exon snap 100 kb. Twenty
independent null replicates calibrate the maximum-statistic threshold. Eight new
null, eight single-sweep, and eight two-sweep replicates evaluate it. Complete
simulated segments are independent units; neither windows nor sites are randomly
split between calibration and evaluation. RAiSD and the composite/gap/linkage
ensemble provide different summary-statistic checks, although both use the same
simulated genotypes. The
[RAiSD study](https://www.nature.com/articles/s42003-018-0085-8) likewise separates
null calibration from power evaluation and shows why demographic and BGS
backgrounds matter.

At SLiM scaling factor 10, the null threshold is 0.7442. There is one false
positive segment among eight null evaluations. No selected sites are detected
in either eight single-sweep or eight two-sweep evaluations (0/24 selected sites).
Unconditional single-sweep localization is better than detection power: the
three-estimate ensemble has median error 40.0 kb and 90th-percentile error
79.44 kb; RAiSD's global peak has median error 36.02 kb and 90th-percentile error
76.43 kb. Each lies within 80 kb in 7/8 replicates. A reliable location conditional
on knowing a sweep exists does not imply a reliable detector of an unknown count.

The [prespecified design](../experiments/scientific_2026_11_20/sweeps/design.json)
and [benchmark record](../experiments/scientific_2026_11_20/sweeps/background_benchmark.json)
contain all calibration maxima and individual evaluation outcomes. The [exact-scale design](../experiments/scientific_2026_11_20/sweeps/design_exact.json)
and [exact-scale benchmark](../experiments/scientific_2026_11_20/sweeps/background_benchmark_exact.json)
repeat the same 44-replicate design at Q=1 with new independent seeds and the
same detector settings. No evaluated threshold or width was reused as a tuned
hyperparameter: each scale has its own independent null threshold under the
unchanged order-statistic rule.

| Independent evaluation at Q=1 | Result |
|---|---|
| Null segments | 0/8 have a detected region; the two-sided exact 95% interval for the false-positive probability is 0–36.9%. |
| One sweep | 7/8 sites covered (exact 95% interval 47.3–99.7%), but only 7/13 reported regions contain a true site (53.8% region precision). |
| Two sweeps | 8/16 sites covered (50% recall); 8/11 reported regions contain a true site (72.7% region precision). |
| Three-estimate single-sweep localization | Median error 37.89 kb, 90th-percentile error 170.24 kb; only 5/8 within 80 kb. |
| RAiSD global-peak single-sweep localization | Median error 51.23 kb, 90th-percentile error 99.16 kb; 6/8 within 80 kb. |

The Q=10 null false-positive estimate of 1/8 likewise has a very wide exact
95% interval, 0.3–52.7%, and does not establish a violation of nominal 5%
control. These [uncertainty calculations](../experiments/scientific_2026_11_20/sweeps/uncertainty.json)
use independent segments; multiple sites or regions within one segment are
correlated and must not be counted as independent binomial trials.

The contrast with Q=10 shows that the approximation can change practical
conclusions under this BGS/selection/history configuration. Separate seeds mean
this is an independent distribution comparison, not a paired equivalence test.
The exact-scale result is materially better for detection but still misses half
the selected sites in two-sweep cases. Neither scale supplies a justified
Challenge 7/8 freeze.

## Remaining scientific gates

None of Challenges 5–8 is certified for final submission by these checks.
Challenge 5 needs fresh full-chromosome recovery, callable-mask and map-domain
stress, and evidence beyond the small consensus-selected subset. Challenge 6
needs a better independently developed inference rule, fresh evaluation of weak
and incomplete standing-variation sweeps, and evidence for other soft-sweep
mechanisms if those are in scope. The current standing-variation simulator is
not a recurrent-mutation soft-sweep generator.

Challenge 7 needs recovery of its exact ten-estimate workflow under a varied BGS
and bottleneck/growth background. Challenge 8 needs adequate detection power,
localization and count/overlap behavior, whole-chromosome null maxima, sweep
separations and ages varied across their domain, independent parameter-shift
strata, longer BGS burn-in, and a properly validated scaling approximation.
The new benchmark uses a single 1 Mb region, fixed nuisance parameters and sweep
locations, few evaluation replicates, and burn-in 1. It cannot establish genome
wide error control or steady-state BGS adequacy. The exact-scale follow-up checks
one scaling concern without removing these limitations. No threshold was tuned
to the evaluation results or any public challenge answer.

## Reproduction and execution

The runner is `python -m scripts.validate_sweep_background_smoke`; its output
contains the scenario, seed and report identity. The scorer is
`python -m scripts.score_sweep_background_validation --input-dir INPUT
--coding-bed EXONS --design DESIGN.json --out RESULTS.json`. The hard/soft
reanalysis uses `python -m scripts.recompute_sweep_recovery --record OLD.json
--trees SIM.trees --scan SCAN.csv --out NEW.json`, followed by the existing
recovery summarizers. Cluster batch files, software identities, source hashes
and job outcomes are retained in
[verification](../experiments/scientific_2026_11_20/sweeps/verification.json).
The executed generator and localization helper differ from the repository only
in formatting, with identical Python AST hashes. All 28 focused sweep tests and
critical Ruff checks pass. These are software checks, separate from the
scientific failures and limitations above.
