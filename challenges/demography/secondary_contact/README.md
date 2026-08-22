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

The exact full-resolution SFS fits are stable to polarization. Folded and
unfolded estimates respectively give mainland sizes of `18,725` and `19,041`,
island sizes of `51,450` and `50,733`, ancestral sizes of `16,578` and `16,593`,
split times of `4,638` and `4,607` generations, migration-start times of `1,183`
and `1,193` generations, and migration rates of `2.45e-4` and `2.40e-4`.
Projections to 40×30 and 20×15 chromosomes land in the same local basin and are
used as sensitivity checks rather than independent estimates.

In addition to learned genealogy features, a continuous-time three-state
coalescent reference now propagates mainland/mainland, cross-population, and
island/island lineage-pair states through contact, isolation, and the ancestral
epoch. It fits span-weighted POLEGON interval counts directly and provides an
interpretable guard for the learned posterior. On the first known-parameter
calibration, POLEGON improved the global hazard RMSE from `0.0444` to `0.0400`
and raised all pair-class correlations above `0.994`, so the redated branch is
promoted for chromosome-wide extraction.

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
