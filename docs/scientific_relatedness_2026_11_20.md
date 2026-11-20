# Scientific review of relatedness challenges 9–12

These results validate components of the estimators, not perfect inference or a
complete final submission. All experiments below use new simulation seeds. No
public challenge answers, final inputs, or submissions were used. Models were
held fixed except for a separately identified Task 9 candidate refitted only to
remove retired labels from its training contract.

## Design and independence

The initial Betty array, job 10849217, generated independent Ped-sim datasets for
Task 9 with sex-specific maps and both interference and Poisson crossover
models, human family relationships with Poisson crossovers, and dog relationships
with the average dog map and Poisson crossovers. Three independent seeds were
used per setting. Job 10849361 adds three human datasets with the matched
interference model. The definitions are the existing independent-pedigree
calibration templates, not whole 229- or 317-person challenge families. New
seeds test inheritance variability and crossover-model transfer; they do not
establish transfer to new pedigree topology or founder population history.

The canonical results are the JSON files under
[`experiments/scientific_2026_11_20/relatedness`](../experiments/scientific_2026_11_20/relatedness/).
The job scripts, model/input SHA-256 hashes, and fixed seeds make the evaluation
reproducible. The scorer performs no fitting or model selection. Its interval
resamples entire independent simulated pedigrees, keeping multiple pairs from
the same pedigree together. Background intensity settings and map perturbations
on the same pedigree are repeated measurements, not additional independent
samples. The 95% intervals therefore describe inheritance sampling for these
fixed template mixtures, not uncertainty over biological model assumptions.

The final evidence contains 17 independent simulation seeds and 25 component
result files. Scoring the same dataset with two models does not create another
independent replicate. `summarize.py` beside those files rebuilds the aggregate;
`verification.json` and `slurm_accounting.psv` record source hashes and job status.

The model inputs are exact Ped-sim IBD segments. Hap-IBD detection, phasing,
genotype errors, correlated founder background, dog inbreeding, HBD mask
fragmentation, family clustering and final pedigree consensus are not executed
by this experiment. The Ped-sim documentation explicitly distinguishes direct
IBD simulation from genetic-data simulation requiring founder genotypes. Crossover
interference and sex-specific maps affect IBD sharing distributions, which is
why both crossover models are tested. These are established properties of the
[Ped-sim method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007979)
and its [official implementation](https://github.com/williamslab/ped-sim).

Hap-IBD uses phased genotypes, genetic-length seed thresholds and gap extensions;
its detection errors depend on marker informativeness and genotype error.
Exact-segment results cannot certify the downstream distribution of detected
segments. See the [hap-IBD paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7118582/)
and [official parameters](https://github.com/browning-lab/hap-ibd).

## Measured component performance

All values below condition on at least one true IBD segment. Exact accuracy for
the frozen Task 9 model counts zero accepted-probability failures as wrong.
Numbers are the mean across independent seeds; map/background variations within
a seed are not counted again. ECE is the absolute reliability error in ten
equal-width confidence bins, not a guarantee of posterior calibration.

| Task and condition | Seeds | Exact accuracy | Meioses accuracy | ECE |
|---|---:|---:|---:|---:|
| 9, frozen; interference | 3 | 88.87% | 89.20% | 1.12 pp |
| 9, frozen; Poisson shift | 3 | 87.21% | 87.83% | 1.66 pp |
| 9, ten-label candidate; interference diagnostic sets | 3 | 88.92% | 89.25% | 1.20 pp |
| 9, ten-label candidate; Poisson diagnostic sets | 3 | 87.26% | 87.86% | 1.78 pp |
| 9, ten-label candidate; untouched interference confirmation | 1 | 89.10% | 89.41% | 0.88 pp |
| 9, ten-label candidate; untouched Poisson confirmation | 1 | 87.20% | 87.71% | 1.97 pp |
| 10, matched interference | 3 | 47.50% | 69.28% | 2.52 pp |
| 10, Poisson shift | 3 | 44.78% | 65.79% | 5.28 pp |
| 11, dog average map / Poisson | 3 | 63.17% | 69.94% | 2.16 pp |
| 12, matched interference / 0.6× empirical contamination | 3 | 38.49% | 55.26% | 2.26 pp |
| 12, matched interference / 1.4× empirical contamination | 3 | 36.36% | 51.96% | 3.02 pp |

These point averages conceal important strata. On the first Task 9 matched seed,
balanced accuracy is 83.35%, versus 88.97% pair-weighted accuracy; fifth-cousin
recall is only 37.42% among pairs sharing any IBD. Roughly 1,000 of 11,000
generated accepted-label pairs are absent because they share no IBD. Human
templates retain 6,743–6,817 of 7,000 pairs in these runs. A workflow requiring
every genotyped pair must also validate the omitted zero-IBD cases.

The historical forest rejects three unperturbed Poisson pairs across the three
new datasets; any one such pair aborts the strict production predictor for that
input. Map-length stresses reveal additional failures. The ten-label candidate
uses 60,266 original calibration pairs. Its two untouched confirmation seeds
(2026112307 and 2026112308) were generated only after the fitted model was fixed.
They show no zero accepted-probability failures under any of the three map-length
conditions. Confirmation exact accuracy is 89.10% versus 89.10% for the frozen
forest under interference, and 87.20% versus 87.06% under Poisson. This supports
the probability-support repair, not a substantial accuracy improvement.
The candidate also has no such failures across the six original diagnostic
datasets and their map-length stresses. Those six datasets first exposed the
problem and are reported separately from the untouched confirmation.

The ±10% map stress scales segment and sex-map genetic lengths while keeping
the normalization genome length fixed. It probes a measurement-scale mismatch;
it is not a newly generated recombination map or a frequency-weighted model for
map uncertainty.

## New failure modes and corrections

**Task 9's retired-class model can reject valid current-domain inputs.** The
historical forest contains twelve classes. Renormalizing it over the current ten
accepted classes is undefined when all its trees assign the retired classes.
Fresh Poisson runs exposed that failure; the strict predictor correctly rejects
such rows. The scientific scorer reports their count, conditional performance
on scorable rows, and accuracy counting failures as wrong. It never inserts a
fallback relationship or silently drops a failed row from the headline result.

`train_task9_current_candidate.py` clones the frozen forest's hyperparameters,
refits on original calibration replicates 0–5 with only the ten currently
accepted labels, and preserves provenance. It performs no search against the
new evaluations. The candidate is not promoted or used for challenge predictions.
The ten-class refit addresses probability support; it does not establish better
biological transfer. Task 9 additionally has fixed ±10% map-length perturbations,
a mechanistic Mendelian expected-sharing reference, classwise recall, log loss,
Brier score and reliability bins. The expected-sharing comparison groups
indistinguishable expectations and is not advertised as an exact-relation
classifier. Meioses and kinship are not interchangeable: for example, full and
half first cousins have the same number of meioses but different expected sharing.

**Task 12's previous background experiment had an unavailable clean reference.**
It corrected each contaminated evaluation replicate against the very same clean
pairs used to generate it. This is not public-answer leakage, but it supplies a
counterfactual unavailable at prediction time and understates reference
uncertainty. The empirical training pipeline now uses other calibration
replicates as its correction reference, leaving out the current replicate.
Selection and development use calibration-only references; the holdout CLI
requires a distinct explicit reference and rejects content-identical inputs.
The new scientific evaluation uses original calibration replicate 0, independently
of every newly generated evaluation seed. Repeated contamination regimes are
reported as correlated measurements; the positivity gate now counts independent
pedigree replicates. The holdout scorer now also honors the saved classifier
blend weight; the historical candidate weight was 1.0, so that correction does
not change its old score. The freeze writer refuses legacy summaries lacking the new
independent-reference evidence. Historical freeze files are preserved, not
recertified by this review. Retraining under the corrected protocol is still
required before promotion.

The older generic background validator now also requires a separate reference.
Its seeded individual effects now iterate through sorted sample IDs, eliminating
Python hash-order dependence. The training entry points reject duplicate files
across calibration/selection/evaluation, including renamed identical copies.
Hashes do not by themselves prove independent simulation seeds or genealogies;
that still needs recorded provenance. Histogram boosting's internal random-row
early stopping is disabled in the relatedness trainers because it can split
correlated pedigree pairs or augmented copies. External simulation-level
selection remains in place. None of these changes retroactively retrains the
historical model artifacts.

Task 11's historical freeze plot now handles permitted labels absent from a
cohort and derives the family count from its input rather than assuming eight.

## Interpretation and remaining validation

Task 9 is the strongest of these components but remains conditional on observing
IBD and the simulation model. Segment files omit pairs sharing no IBD, especially
for distant cousins. We report expected generated pairs and retained counts so
this conditioning is visible. The ten-label refit now has untouched
confirmation, but still requires a release decision and current-input validation
before use. Current output-format
validity is separate from relationship accuracy.

Tasks 10–12 do not yet satisfy the full validation protocol. Their calibration
uses isolated relationship templates; the complete challenge estimator also
uses all-pair topology, HBD masking or cohort-derived background estimates.
New exact-segment results cannot establish the accuracy of these composed
procedures. Pedigree reconstruction stability across settings is not accuracy,
and averaging its outputs does not supply calibrated posterior probabilities.

Task 10 requires joint large-family simulation through genotype generation,
phasing/IBD detection and the frozen consensus, holding out entire families.
Task 11 additionally needs matched close-inbreeding founder pedigrees, the same
HBD masking and family clustering, and evaluation of missed/merged families.
Task 12 needs independently generated background genealogy, varied family
composition and a separate reference cohort. Its frozen empirical bin profile
and Poisson-added segments are a diagnostic contamination model; they cannot
reproduce cross-pair dependence induced by shared founders. Probability
calibration and uncertainty must be assessed after these full pipelines, not
borrowed from the pairwise forest or boosting classifier.
