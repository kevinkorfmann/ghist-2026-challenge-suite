# Task 9 methods — final run, November 20, 2026

The frozen ten-label release was run on Synapse `syn75202631`, current version 3, on November 20, 2026, in Slurm job 10856122 (23 seconds). All 997 observed pairs from 39,916 segments are represented once. An independent coordinator audit confirms all ten run-file hashes, complete pair coverage and valid relationship/meiosis assignments. Final accuracy is unknown.

Input SHA-256: `ff7031bc00a96f7d38131e010f1168f95fb4a8efad5849c904b2d7e77e4de394`. Prediction SHA-256: `d2810e81b74591c3d10c5305643c3ed4494ab416e44e4f98c37ff6df92b16c6d`. Uploaded as [syn78778815](https://www.synapse.org/Synapse:syn78778815), version 1, and submitted as an individual to GHIST 2026 Relatedness 1 Final. Synapse displayed its successful submission receipt. Scoring is pending.

We classify relationships from exact pairwise Ped-sim IBD1 and IBD2 segments. The
operational input domain is the set of distinct unordered individual pairs
present in the segment CSV. The supplied specification, as recorded in the
[readiness report](task9_readiness_2026_11_20.md), does not identify a separate
roster of zero-IBD pairs. Consequently, coverage of the observed pair domain is
verified, and no unobserved pair or identity is invented. Accuracy estimates are
conditional on a pair having observed IBD. This domain interpretation remains a
limitation if the challenge expects predictions for pairs absent from its CSV.

Each pair is represented by 218 fixed features. These include IBD1/IBD2 segment
counts and total lengths; segment-length moments, quantiles, histograms and
ordered lengths; chromosome-level sharing summaries; IBD2 proportions; and
male/female map-derived segment lengths and their differences. Physical segment
endpoints are interpolated on the frozen sex-specific map. Relatedness is
normalized by a genome length of 3,350.933402 cM. The release manifest fixes the
feature order, map checksum and feature/label source checksums.

We fitted a random forest to 60,266 pairs from the six original calibration
simulation replicates. Here, “calibration” denotes the training simulations;
classifier probabilities were not subsequently calibrated. The refit removed
retired training labels while retaining the historical forest settings: 600
trees, Gini splitting, bootstrap sampling, balanced-subsample class weights,
`max_features=0.8`, minimum leaf size 2, and random seed 90926. Other settings are
recorded in the [training report](../experiments/scientific_2026_11_20/relatedness/task9_current_training.json).
No public answers or new confirmation outcomes were used to select the refit's
hyperparameters. The refit repaired accepted-label probability support; it is
not evidence of a substantial improvement over the historical forest.

The ten accepted output labels and deterministic meiosis assignments in this
release are:

| Relationship label | Meioses count |
|---|---:|
| `siblings` | 2 |
| `grandparent_grandchild` | 2 |
| `avuncular` | 3 |
| `double_cousin` | 4 |
| `full_1st_cousin` | 4 |
| `half_1st_cousin` | 4 |
| `full_2nd_cousin` | 6 |
| `full_3rd_cousin` | 8 |
| `full_4th_cousin` | 10 |
| `full_5th_cousin` | 12 |

The fixed model was evaluated on two subsequently generated, untouched
confirmation datasets. The matched interference scenario used seed 2026112307;
the Poisson-crossover shift used seed 2026112308. These simulations were distinct
from the six calibration files by content checksum. Exact accuracy is the
fraction of correctly classified relationship labels, balanced accuracy is the
mean label recall, and meioses accuracy compares the mapped counts.

| Untouched confirmation | Pairs | Independent pedigrees | Exact accuracy | Balanced accuracy | Meioses accuracy |
|---|---:|---:|---:|---:|---:|
| Matched interference | 10,038 | 9,038 | 89.10% | 83.35% | 89.41% |
| Poisson crossover shift | 9,985 | 8,985 | 87.20% | 81.48% | 87.71% |

Pedigree-bootstrap 95% intervals for exact accuracy were 88.45–89.70% and
86.53–87.87%, respectively. The preserved minimum release gates required at least
10,000 matched pairs, 80% exact accuracy, 75% balanced accuracy and 84% meioses
accuracy; both confirmations passed the numerical accuracy gates. Additional
map-length stresses retained valid probability support. These gates concern
the tested simulations and do not establish accuracy on current final data.

The actual production predictor and a separately invoked submission-validator
CLI were exercised on anonymized segment CSVs with independently recorded pair
lists. Equivalent row/column orderings, reversed pair order and UTF-8 BOM inputs
gave the same pair-to-label predictions and reproduced the fixed confirmation
scores. Eight production executions, including the explicitly versioned release
path, passed across 20,023 distinct confirmation pairs. Input checks reject
duplicate segment records (including reversal of the same pair), HBD,
malformed/nonfinite or negatively valued coordinates, reversed intervals,
inconsistent genetic lengths, missing map chromosomes, empty/nonfinite or
nonmonotonic maps, and pair-coverage discrepancies. These checks do not establish
biological validity: the frozen feature code permits fractional physical
coordinates, zero-length segments and overlapping intervals. Physical endpoints
outside a map's represented range are clamped to its endpoints by NumPy
interpolation; map-relative bounds are not rejected. Output consists of one row per
observed unordered pair with columns `indID1,indID2,meioses_count,relation`.

The model remains limited by the simulated pedigree templates and the tested
crossover/map scenarios. Distant-cousin distinctions are difficult: matched
confirmation recall was 56.66% for fourth cousins and 33.87% for fifth cousins.
Forest class probabilities are not calibrated posteriors. This validation uses
exact segments and does not certify genotype processing, phasing, detected-IBD
errors, HBD, background ancestral sharing, or complete-pedigree reconstruction.

The inference handoff package pins the complete prediction/validation source,
release manifest, evidence files and runtime requirements. The actual final-input identity, output identity and coverage audit are recorded above and in the retained run handoff.
The package does not contain the calibration/confirmation simulations or their
generators and therefore is not a self-contained retraining package. Its evidence
checks establish consistency with the recorded validation; they do not rerun that
validation or authenticate externally supplied Synapse provenance.
