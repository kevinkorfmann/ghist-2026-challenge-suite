# Independent review of the complete-pipeline validation wrappers

The complete Task 10 and Task 11 simulations exercise generated genotypes,
preparation, real hap-IBD detection, feature extraction, clustering, frozen pair
models and pedigree consensus. Their truth is used for scoring, not inference.
The missingness experiment deliberately stops at pair classification and does
not claim to rerun the already failing consensus. Supplied family count, simple
founder models, true phase and restricted relationship classes remain explicit
limitations. None of these studies qualifies a submission.

An independent read-only recalculation reproduced Task 10's 58 correct labels
among all 153 pairs, including reconstruction of truth from the retained full
pedigree and identifier mapping. Both Task 11 cohorts retain all 48 true pairs
in their denominators: consensus has 7 and 9 correct labels, respectively. All
three frozen-forest thresholds, family partitions and eight missingness policy
scores also reproduce their stored results. The eight policies reuse the same
96 distinct pairs from two cohorts; they are not new biological replicates.
All 42 locally retained files covered by Task 11's recorded output manifests
match their hashes. Fourteen larger files are remote only and were not rehashed
by this local review.

The review found a reproducible identity gap in the missingness wrapper. It
checked baseline genotype bytes and seed but did not compare `truth.json` or
`protocol.json` with the hashes recorded by the baseline run. Replacing the
scoring relationship while preserving the seed and VCF was accepted and scored
under the substituted truth. This was reproduced using temporary synthetic
fixtures and mocked expensive inference stages, not by altering experiment
records. No evidence was found that the retained studies suffered such a
replacement.

New wrapper executions now require the baseline truth, protocol and VCF hashes
to match the original result, and separately verify the result's protocol
linkage. Task 11 and missingness recheck recorded sources, models, maps,
detector artifacts and generated scoring/genotype inputs before accepting
completed inference. Task 10 already rechecked its source and external inputs.
Task 10 and missingness now record all retained generated output hashes,
including predictions, prepared genotypes and truth where generated, alongside
the exact protocol hash. All three use the explicit
`relatedness-pipeline-v2` provenance schema. A failed check remains a failure;
finite predictions or completed external commands cannot override it.

The original experimental records, scores, source snapshots and verification
files remain unchanged. Task 10 and missingness originally lacked the new
output manifests; Task 11 and missingness did not have the new end-of-run
input checks. Those historical runs must not be described as having passed the
new safeguards. The current regressions validate the safeguards themselves,
not improved biological performance or a new complete-pipeline simulation.

Validation: 26 focused tests pass across
`test_task10_complete_pipeline.py`, `test_task11_complete_pipeline.py` and
`test_relatedness_missingness.py`; Ruff passes for those tests, the three
wrappers and `src/ghist_relatedness/provenance.py`. Added regressions cover
same-seed truth substitution, missing truth digests, changed protocols and
models, post-generation truth changes, and exact candidate/output bindings.
An expanded compatibility check including genotype preparation, haplotype
conversion, Task 11 consensus and its freeze guard also passes all 180 tests.
No challenge answers or final challenge data were read, and no submissions
were performed.
