# Secondary contact

Infer five physical parameters for a mainland–island split followed by delayed,
continuous migration:

1. mainland effective population size,
2. island effective population size,
3. generations since the split,
4. generations since migration began, and
5. per-generation migration probability.

The testing VCF contains 43 mainland and 24 island diploid samples across one
100 Mb contig. The registered mutation and recombination rates are `8.4e-8` and
`5.39e-8` per base per generation. The exact template and score convention are
frozen before fitting. Final files are downloaded but sealed outside this
challenge directory and are not inspected during method selection.

The exact testing audit retained 4,812,034 biallelic records and filtered 83,302
non-biallelic records. Mainland and island diversity per base are `0.0060065`
and `0.0061737`; cross-population divergence is `0.0063144`, with Hudson's
`FST = 0.03552`. All registered ancestral-allele values equal REF, so folded and
unfolded fits are both retained as a polarization-sensitivity check rather than
assuming the annotation adds independent orientation information.

## Method stack

The mechanistic references are folded and unfolded two-dimensional SFS fits in
moments and dadi, global-plus-local GADMA search, and joint SFS+LD inference in
momentsLD. The learned branch combines a two-dimensional SFS encoder,
multiresolution chromosome tokens, within-population LD, private and shared
allele spectra, divergence summaries, and cross-population coalescence hazards.
Score-aligned stacking is eligible only after parameter-point-isolated selection
and untouched whole-chromosome evaluation.

The genealogy branch samples posterior tree sequences with
[SINGER](https://github.com/popgenmethods/SINGER), re-estimates branch lengths
with [POLEGON](https://github.com/YunDeng98/POLEGON), and converts three pair
classes—mainland/mainland, island/island, and mainland/island—into conditional
coalescence hazards on a fixed log-time grid. Posterior draws are averaged only
after calibration against known simulated trees; topology-only SINGER,
SINGER+POLEGON, and truth-tree hazards remain separate ablations.

Two independent full-resolution folded SFS searches give mainland sizes of
`18,725` and `19,041`, island sizes of `51,450` and `50,733`, ancestral sizes of
`16,578` and `16,593`, split times of `4,638` and `4,607` generations,
migration-start times of `1,183` and `1,193` generations, and migration rates of
`2.45e-4` and `2.40e-4`. Projections to 40×30 and 20×15 chromosomes land in the
same local basin. A mode-label audit found that an earlier shell condition also
folded the run labeled as unfolded; corrected polarized fits are therefore kept
separate and are not claimed until their new jobs finish.

The pinned historical generator exposes deterministic size epochs that a simple
five-parameter SFS omits: three ancestral changes just before the split, two
mainland changes after it, and two island changes around contact. A second
mechanistic model now integrates those epochs exactly, parameterizing only the
five reported targets plus the ancestral scale inferred from theta. Folded and
correctly polarized searches at three projections are compared to the simple
model before either can be promoted.

The same exact epoch convention is also used by a dedicated score-aligned SFS
training campaign. This campaign is isolated from the earlier simple-model
surrogate so that its untouched recovery score measures the correction rather
than averaging incompatible generators. The earlier branch remains a negative
control until the generator-matched selection and guard results are complete.

A separate full-chromosome parametric bootstrap matches the observed finite-site
mutation process with discrete genomic positions and JC69 mutations before
filtering multiallelic records. Its smoke replicate produced a `1.67%`
multiallelic fraction, compared with `1.70%` in the testing VCF. Thirty-two
independent 100 Mb replicates are used to estimate local finite-genome bias and
sampling intervals for the five reported parameters.

On an independent known-parameter 10 Mb recovery case, the exact hidden-epoch
model achieved competition RRMSE `0.309`, versus `3.195` for the simple model.
The matched fit reduced the error by more than tenfold, while retaining a
contact-time/migration tradeoff that the full-chromosome bootstrap must quantify
before any testing estimate is frozen.

In addition to learned genealogy features, a continuous-time three-state
coalescent reference now propagates mainland/mainland, cross-population, and
island/island lineage-pair states through contact, isolation, and the ancestral
epoch. It fits span-weighted POLEGON interval counts directly and provides an
interpretable guard for the learned posterior. On the first known-parameter
calibration, POLEGON improved the global hazard RMSE from `0.0444` to `0.0400`
and raised all pair-class correlations above `0.994`, so the redated branch is
promoted for chromosome-wide extraction.

The inference-error calibration is not based on one favorable example. A
separate grouped campaign simulates 768 parameter points with three independent
100 kb replicates each, runs the same SINGER and POLEGON settings used for the
testing chromosome, and trains only on the resulting inferred hazards. Core,
guard-shell, and broad-prior points remain labeled through untouched evaluation.
The observed chromosome is divided into 1,000 non-overlapping tiles and is
aggregated through fifty 2 Mb blocks, preserving block uncertainty and preventing
any single inferred genealogy from dominating the result.

The pinned historical GHIST generator is used only as a convention and
misspecification audit. Its secondary-contact data included unreported ancestral
and population-specific size epochs around the split and contact. Validation
therefore randomizes comparable hidden epochs, directional migration, local
recombination variation, mutation-rate error, missing calls, ancestral-state
error, and smooth linked-selection suppression. Half of the training mass stays
strictly nominal for precision.

The full frozen protocol is in
[`configs/secondary_contact_experiments.yaml`](../../../configs/secondary_contact_experiments.yaml).
No prediction is uploaded or submitted automatically.
