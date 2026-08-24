# Method comparison and extension plan

The project uses mechanistically distinct references and learned estimators. A
method advances only when it improves held-out accuracy without degrading
calibration or robustness under misspecification.

## Demographic history

| Family | Representative method | Signal | Role here |
|---|---|---|---|
| Diffusion/composite likelihood | [dadi](https://doi.org/10.1371/journal.pgen.1000695) | one- or multi-population SFS | independent diffusion cross-check |
| Moment equations | [moments](https://doi.org/10.1534/genetics.117.200493) | SFS moments | primary mechanistic reference |
| Linkage moments | [momentsLD](https://doi.org/10.1093/molbev/msz237) | LD statistics through distance | timing-sensitive reference |
| Simulation composite likelihood | [fastsimcoal2](https://doi.org/10.1371/journal.pgen.1003905) | simulated expected SFS | flexible complex-model reference |
| Sequence-position embedding | [SPIDNA/dnadna](https://doi.org/10.1093/bioinformatics/btac765) | genotype matrix and SNP positions | raw-input comparison |
| Posterior estimation | [Smith et al. 2026](https://doi.org/10.1093/genetics/iyag107) | raw genotypes or SFS+LD | calibrated posterior target |
| Genealogical hazards | [Deng et al. 2026](https://doi.org/10.64898/2026.08.11.744210) | interval-conditional coalescence | temporally localized extension |
| Coalescence translation | [CXT](https://github.com/kr-colab/cxt) | genotype-conditioned TMRCA posterior | independent genealogy channel |
| State-space translation | [FastCXT](https://github.com/kevinkorfmann/fastcxt) | bidirectional Mamba log-TMRCA mean and variance | fast sample-size-flexible genealogy channel |

### Extension adopted here

For each finite interval, the genealogical feature is

`P(T in interval | T has not occurred before interval)`.

This removes the temporal-covariance leakage that arises when sites or branches
are merely filtered by mutation age. For bottleneck inference, the resulting
hazard curve should expose the population-size transition directly. For secondary
contact it will be computed within and between populations, where the contrast is
expected to localize split and renewed migration. Genealogical channels are
randomly masked during training and benchmarked under inferred-tree error.

CXT integration uses the pinned upstream revision
`294b0c27226accbac830c9c922de6b2f92f83fe8`. Posterior log-TMRCA samples are
converted to the same conditional hazard bins used for simulated tree sequences,
block by block. This permits a direct comparison among CXT, tsinfer/tsdate, and
mutation/LD-only predictions, while retaining CXT posterior uncertainty instead
of collapsing it to one genealogy.

FastCXT integration uses pinned revision
`7fb89f451cca734e03bb974bca7faffe9a94e194`. For every pair/window Gaussian in
log time, its probability mass in each finite interval is evaluated by the normal
CDF. These masses are averaged and divided by the posterior mass surviving to the
interval start. This is an exact posterior hazard under the model and avoids
Monte Carlo noise. One-megabase translation blocks may be combined into the
two-megabase feature blocks through explicit genomic coordinates. FastCXT is
trained on the challenge prior and nuisance domain before comparison; an
untrained upstream model is never treated as evidence.

The production posterior combines:

- folded SFS and logarithmic LD decay over the full sequence;
- multi-scale window summaries;
- optional conditional coalescence hazards;
- a permutation-aware context encoder;
- a flexible conditional density head and independent simulation-reference trees.

The architecture screen also contains a direct FastCXT-inspired candidate. It
replaces only the global transformer with stacked bidirectional Mamba-2 blocks;
the local dilated convolutions, inputs, mixture-density head, split, training
budget, and promotion metric remain fixed. This isolates the value of state-space
context instead of confounding it with a different feature set or validation
protocol.

## Sweep detection

| Family | Representative method | Role here |
|---|---|---|
| Composite likelihood | SweepFinder2 | likelihood reference |
| SFS/LD scan | RAiSD and OmegaPlus | fast references |
| Spatial summary classifier | S/HIC and diploS/HIC | demographic-robust learned reference |
| Wavelet summaries | SURFDAWave | multi-scale reference |
| Flexible sweep classifier | [Flex-sweep](https://doi.org/10.1093/molbev/msad139) | broad sweep-state reference |
| Raw sequence classifier | [ASDEC](https://doi.org/10.1093/bioinformatics/btad246) | raw-input comparison |

The extension is a chromosome-context segmentation model with accessibility,
coding status, recombination rate, classical scan scores, raw genotype tokens, and
genealogical imbalance. It predicts calibrated per-base sweep probability and an
interval boundary distribution rather than independent window classes.

## Relatedness

Mechanistic references include KING/PC-Relate kinship, hap-IBD and refined IBD,
IBDkin, ERSA-style likelihoods, and pedigree-graph solvers. Learned models operate
on pairwise IBD length spectra, IBD1/IBD2/HBD totals, chromosome counts, sex-map
features, and consistency messages over the full family graph. Validation is held
out by pedigree or family; pairwise random splits are prohibited.

## Decision criteria

1. Freeze priors and validation partitions before viewing public answers.
2. Compare point error, interval coverage, calibration, runtime, and failure strata.
3. Require superiority on independent simulations and at least parity with every
   applicable reference method.
4. Retain ensembles only when disagreement is calibrated and predictive of error.
5. Record all negative results; complexity is kept only when it earns its cost.

## Task 12 decision

Task 12 uses a 5 cM corrected-domain histogram-gradient model. Background rates
are estimated from unlabeled cohort-versus-reference length-bin differences, then
varied across Ped-sim calibration without using testing labels. Against the prior
20 cM model, exact accuracy increased from 34.85% to 37.53% and meiosis-degree
accuracy from 51.91% to 53.24% over 20 untouched empirical-background comparisons.
Every paired comparison improved. Degree-guarded pedigree topology refined 36 of
50,086 frozen testing pairs.

## Task 10 decision

Task 10 uses a 75:25 blend of a shallow histogram-gradient full-shape model and
the prior calibrated forest. The blend improved all eight untouched pedigree
replicates: exact accuracy rose from 46.51% to 47.20% (+0.70 points), while
meiosis-degree accuracy rose from 68.61% to 69.08%. Two independent pedigree
reconstructions may change a call only when they agree exactly and preserve the
pair model's meiosis degree; this refined 198 of 26,106 frozen testing pairs.

## Task 1 decision

Task 1 is frozen locally without submission. Exact folded-SFS moments is primary
at 4,587.06 generations and a 0.5224206 post-decline fraction. It had the lowest
natural-scale medians on the untouched nominal chromosome cells. Independent
dadi and joint momentsLD fits were within 1.3% and 1.6% of its observed decline
time, respectively, and within 0.0018 in fraction.

The five-seed context ensemble, the FastCXT-enriched branch, affine SFS
calibration, and a selection-chosen learned–mechanistic blend were all retained as
diagnostics but rejected for primary use. Each failed either exchangeable input
translation or the frozen natural-error guard. The complete decision record is
[`experiments/bottleneck/task1_freeze.json`](../experiments/bottleneck/task1_freeze.json).

The final precision pass deliberately covered the full held-out prior rather than
only the observed neighborhood. It used 900 complete 100 Mb simulations at 150
parameter points, with all six replicates of each point confined to one split.
The score-aligned two-layer calibrator passed the untouched RRMSE gate, reducing
broad-prior RRMSE from 5.789 to 0.898 and focused RRMSE from 0.506 to 0.206. Its
observed prediction is 4,338.37 and 0.49587. Because it increased typical focused
natural errors and point-mean bias, it is a held-out extreme-value contingency;
the mutually agreeing mechanistic fits remain primary for the present VCF.
