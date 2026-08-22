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

## Method stack

The mechanistic references are folded and unfolded two-dimensional SFS fits in
moments and dadi, global-plus-local GADMA search, and joint SFS+LD inference in
momentsLD. The learned branch combines a two-dimensional SFS encoder,
multiresolution chromosome tokens, within-population LD, private and shared
allele spectra, divergence summaries, and cross-population coalescence hazards.
Score-aligned stacking is eligible only after parameter-point-isolated selection
and untouched whole-chromosome evaluation.

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
