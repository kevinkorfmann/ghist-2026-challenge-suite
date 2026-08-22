# Method comparison and extension plan

The project uses mechanistically distinct references and learned estimators. A
method advances only when it improves held-out accuracy without degrading
calibration or robustness under misspecification.

## Demographic history

| Family | Representative method | Signal | Role here |
|---|---|---|---|
| Diffusion/composite likelihood | [dadi](https://doi.org/10.1371/journal.pgen.1000695) | one- or multi-population SFS | primary mechanistic reference |
| Moment equations | [moments](https://doi.org/10.1534/genetics.117.200493) | SFS moments | fast likelihood reference |
| Linkage moments | [momentsLD](https://doi.org/10.1093/molbev/msz237) | LD statistics through distance | timing-sensitive reference |
| Simulation composite likelihood | [fastsimcoal2](https://doi.org/10.1371/journal.pgen.1003905) | simulated expected SFS | flexible complex-model reference |
| Sequence-position embedding | [SPIDNA/dnadna](https://doi.org/10.1093/bioinformatics/btac765) | genotype matrix and SNP positions | raw-input comparison |
| Posterior estimation | [Smith et al. 2026](https://doi.org/10.1093/genetics/iyag107) | raw genotypes or SFS+LD | calibrated posterior target |
| Genealogical hazards | [Deng et al. 2026](https://doi.org/10.64898/2026.08.11.744210) | interval-conditional coalescence | temporally localized extension |
| Coalescence translation | [CXT](https://github.com/kr-colab/cxt) | genotype-conditioned TMRCA posterior | independent genealogy channel |

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

The production posterior combines:

- folded SFS and logarithmic LD decay over the full sequence;
- multi-scale window summaries;
- optional conditional coalescence hazards;
- a permutation-aware context encoder;
- a flexible conditional density head and independent simulation-reference trees.

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
