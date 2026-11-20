# Haplotype-to-diploid IBD conversion repair, November 20

The shared converter now preserves haplotype identities and labels an interval
IBD2 only when the active matches contain `1→1` together with `2→2`, or `1→2`
together with `2→1`. Two overlapping records are insufficient when they share
either individual's same haplotype. Exact duplicates and overlapping copies of
one edge retain that edge until its last record ends, without adding sharing.
Sorting a sample pair also swaps the corresponding haplotype identities.

The implementation is in
[`reconstruct_bonsai_pedigree.py`](../scripts/reconstruct_bonsai_pedigree.py).
`read_haplotype_segments` now returns explicit `(start, end, haplotype1,
haplotype2)` records. `coverage_intervals` deliberately rejects the old ambiguous
two-coordinate inputs. The three production callers inspected—Task 11 observed
forest feature extraction, Task 11 Bonsai reconstruction and human branch
reconstruction—all call `convert_segments` and receive the corrected behavior;
none invokes a legacy interval-only path.

Physical endpoints remain unchanged. Genetic lengths remain `map(end) −
map(start)`, with no added base or positive-length overlap at a shared endpoint.
Adjacent pieces with the same sharing state merge. The original species-map
length threshold is applied before the optional dog-to-human length scaling.
Zero-cM map plateaus are valid but do not emit zero-length output segments.

Malformed haplotype identifiers, self pairs, nonpositive or reversed physical
intervals, nonfinite reported lengths, invalid thresholds/scales, empty maps,
duplicate/decreasing map positions, and decreasing/nonfinite map cM now raise
explicit errors. Map positions must be positive integers. Chromosomes missing
from the map and intervals extending beyond its support are rejected instead
of being silently clipped by interpolation. This stricter input contract must
be checked for each future dataset; verification on the synthetic cohorts does
not establish coverage for unseen challenge inputs.

## Verification and preserved evidence

The focused converter/caller suite passes **82 tests**, including every subset
of the four possible haplotype edges, both valid full-sharing matchings,
one-to-many sharing, exact and nested duplicate edges, reversed sample pairs,
touching endpoints, map interpolation, threshold/scaling order, zero-pair sample
retention, malformed records, invalid coordinates and invalid maps. Ruff passes.
The parent integration run reports **909 tests passed** across the workspace.

The bounded Slurm audit in
[`ibd_conversion`](../experiments/complete_pipeline_2026_11_20/ibd_conversion/)
compares the original frozen converter and corrected converter against preserved
Task 11 synthetic cohorts with seeds 2026112711 and 2026112811. It evaluates raw
conversion at 0/10 cM thresholds and the existing HBD-masked forest feature
extraction at 2/5/10 cM HBD thresholds. It verifies the historical map, model and
source hashes, reproduces historical feature vectors and forest probabilities,
and checks every retained historical cohort file's hash before and after. No
simulation, fitting, final-data inference, topology or consensus rerun occurs.
Historical results and snapshots remain evidence, not edited output locations.

The first job, 10850841, stopped before Python because the cluster's profile
quota hook reads unset `TMUX` under shell strict mode. The batch bootstrap was
corrected to source the profile before enabling strict mode; this was not an
estimator failure.

Retry **10850864 completed in 22 seconds**, with 3,234,208 KiB maximum reported
RSS. Both synthetic detector outputs passed the stricter map support checks.
The old feature vectors and forest probabilities reproduced at all three HBD
thresholds, and all watched historical files retained their original hashes.
Both attempts' logs and scheduler accounting are retained in the audit folder.

| Unthresholded conversion measure | Seed 2026112711 | Seed 2026112811 |
|---|---:|---:|
| Pairs whose IBD2 changed | 36/48 | 36/48 |
| Old summed IBD2 across pairs, cM | 32,854.81 | 37,333.96 |
| Corrected summed IBD2 across pairs, cM | 23,931.22 | 24,905.90 |
| Largest IBD2 reduction for one pair, cM | 500.27 | 711.61 |

All 72 affected pairs were offspring of related parents; all 24 pairs with
unrelated parents were unchanged. With no cM threshold, the total physical/map
coverage union was preserved to floating-point tolerance (maximum difference
below 8×10⁻¹³ cM per pair), while overcalled IBD2 became IBD1. At the production
10 cM threshold, 72/96 pairs also changed; filtering newly split or merged state
intervals additionally changes retained coverage, as expected from the existing
threshold semantics.

After the existing HBD masking, **58/96, 68/96 and 72/96** pairs had changed
feature vectors at the 2, 5 and 10 cM HBD thresholds respectively. The frozen
forest retained **96/96 correct sibling labels** at every threshold. Its largest
absolute class-probability change was 0.04426. These are repeated measurements
of the same 96 pairs in 16 families, not independent validation replicates.
The audit did not rerun Bonsai or consensus, so no new final-pipeline accuracy
is asserted.

Detailed values, per-pair differences, source hashes and preservation checks are
in [`audit.json`](../experiments/complete_pipeline_2026_11_20/ibd_conversion/results/audit.json),
[`conversion_pairs.csv`](../experiments/complete_pipeline_2026_11_20/ibd_conversion/results/conversion_pairs.csv)
and [`masked_feature_pairs.csv`](../experiments/complete_pipeline_2026_11_20/ibd_conversion/results/masked_feature_pairs.csv).

## Calibration implications and limits

The changed IBD1/IBD2 calls can change total sharing, segment boundaries and
counts, length quantiles, threshold features and Bonsai likelihoods. Historical
accuracy and calibration evidence using the old adapter does not certify the
corrected pipeline. Direct pedsim IBD1/IBD2 training records do not use this
converter and are not automatically changed by the repair; compatibility of
those training summaries with newly corrected detector-derived observations
still requires validation. The separate Task 9 feature/release implementation
was not modified.

The preserved smoke had 96/96 correct forest sibling calls but 16/96 correct
consensus calls. That remains the historical finding. A conversion difference
does not establish that this adapter caused all or any particular wrong
topology, nor does retaining the forest labels validate the new complete
pipeline. The reused, sibling-only, independent-marker simulations are a
diagnostic comparison, not independent confirmation or calibrated challenge
accuracy. Heterogeneous inbred families, matched founder structure and missing
genotypes remain outstanding. No submission readiness is claimed.
