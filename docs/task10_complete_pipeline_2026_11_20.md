# Task 10 independent complete-family operational validation

This experiment adds jointly generated genotypes, detected IBD, branch discovery,
two frozen pedigree reconstructions and final guarded consensus to the earlier
isolated exact-IBD pair-template evaluations. It does not by itself qualify a
229-person challenge submission. No final challenge inputs, public answers,
model fitting or submissions are used.

## Fixed design and independence

The available Betty assets were audited first. Ped-sim, hap-IBD, Bonsai, the
sex-specific human map, msprime/tskit and the frozen Task 10 models are installed.
The targeted audit found no independent founder genotype panel suitable for this
test; pre-existing relatedness VCFs are challenge-derived and were not used.
Instead, the experiment generates 14 diploid founders independently of all
training and challenge genotypes, with msprime constant population size 10,000,
mutation rate 1.25e-8, and chromosome-specific mean recombination rates. It uses
all 22 autosomes. Common founder variants (minor allele count at least two of
28 haplotypes) are retained at least 500 bp apart. This ascertainment and simple
founder population history are fixed diagnostic assumptions, not a fitted
reproduction of the challenge population.

Ped-sim propagates these founders through one four-generation tree under the
existing sex-specific map and interference parameters. The terminal generation
contains 18 genotyped people in nine sibling pairs: all 153 possible pairs are
scored, consisting of nine sibling, 36 first-cousin and 108 second-cousin pairs.
Truth is derived from the generated full pedigree, not from pair-name templates
or inference output. The seed is 2026112010; Ped-sim uses seed 2026121010.

The generator retains true phase and disables genotype errors and missingness.
The current official challenge contract states that VCFs are phased and
error-free but can contain missing genotypes. True phase is therefore a matched
assumption, not evidence that a separate phasing algorithm must be added. This
assay covers the missingness-free subset; handling and evaluating missing calls
remain separate work. See the [official relatedness task specification](https://www.synapse.org/Synapse:syn74803798/wiki/640950).

Hap-IBD uses the existing settings: seed 1 cM, extension 0.2 cM, output 2 cM,
100 markers and 1,000-bp gap. The classifier remains the frozen 75% clean-spectrum
boosting / 25% forest blend. Branch discovery retains edges with relatedness
from segments at least 20 cM of 0.12 or greater. Both existing Bonsai branch
settings run: maximum add degree / append types 3/1 and 4/2, retaining segments
at least 10 cM. Consensus requires replicated topology and pair-model meiosis
agreement. The real pedigree is never passed to discovery, reconstruction or
consensus.

The job records all source, model, map and binary hashes and the complete design
before genotype generation; it verifies them again after scoring. The source
snapshot is isolated under Betty `artifacts/task10-complete-20261120/source`.
Full synthetic genotypes and detector files remain in the adjacent seed
directory. [The job script](../experiments/complete_pipeline_2026_11_20/task10/run.sbatch)
is versioned locally. The initial scheduled job is **10850569**.

## Operational contract and readiness

The harness requires all 153 unique, non-self pairs and the official three-column
candidate schema. Missing predictions count as errors in the full-family
denominator; absent zero-IBD pairs fail the coverage contract rather than being
dropped from reported accuracy. Unsupported labels, unknown pairs, duplicate
pairs and changed dependencies fail explicitly. Seven targeted local tests
exercise pedigree-derived truth, all-pair coverage, duplicate rejection and
phased/missing-genotype checks.

The operational pass criterion was fixed before running: every stage completes
and produces all 153 valid pairs. Accuracy is descriptive, with no threshold
selected against this family and no refitting. There is one independent family,
not 153 independent replicates; pair-level confidence intervals would misstate
the evidence.

Job **10850569 completed successfully in 8 minutes 24 seconds**, with maximum
reported resident memory 3,748,316 KB. The family VCF contains 2,229,944 markers;
hap-IBD retained 2,090,206 markers and reported 5,240 IBD and 63 HBD segments.
All 153 pairs reached the final candidate, so the operational coverage and schema
contract passed. Relationship accuracy was poor:

| True relationship | Pairs | Correct exact label | Correct meioses |
|---|---:|---:|---:|
| Siblings | 9 | 5 (55.56%) | 5 (55.56%) |
| First cousins | 36 | 16 (44.44%) | 36 (100%) |
| Second cousins | 108 | 37 (34.26%) | 76 (70.37%) |
| **All pairs** | **153** | **58 (37.91%)** | **117 (76.47%)** |

The flat classifier and final consensus produced the same calls. Branch sizes
were 6, 4, 2, 2, 2 and 2. Five calls had replicated topology support; none was
changed by sibling-component propagation. The consensus's reported stability
of 1.0 is uninformative here: it received a single pair-prediction configuration.
It is neither a measured accuracy nor a calibrated posterior probability.

These are measurements of the **original IBD-to-Bonsai converter and original
consensus** in this run's immutable snapshot. The converter discarded haplotype
indices and classified overlapping pair segments as IBD2 based on coverage;
the separately identified converter repair is not evaluated here. Both original
files are archived as `.py.txt` alongside the results, with hashes matching the
preregistered design. The poor flat-classifier result precedes that converter,
so its measured accuracy cannot be attributed solely to the topology adapter.
No model or threshold was changed after seeing these results.

The measured [result JSON](../experiments/complete_pipeline_2026_11_20/task10/seed_2026112010/result.json),
[stage commands and return codes](../experiments/complete_pipeline_2026_11_20/task10/seed_2026112010/commands.json),
[fixed design and hashes](../experiments/complete_pipeline_2026_11_20/task10/seed_2026112010/preregistered_design.json)
and [Slurm accounting](../experiments/complete_pipeline_2026_11_20/task10/slurm_accounting.psv)
are retained locally. The source snapshot and large synthetic files remain on
Betty. The original environment, including an unmodified Bonsai checkout, is
recorded in [environment.json](../experiments/complete_pipeline_2026_11_20/task10/environment.json).

**Readiness effect: execution and coverage are now demonstrated for one complete
small family, but Task 10 remains unqualified for submission.** Remaining scope
includes many independently held-out families, challenge-scale family size,
distant and removed relationships, founder-model transfer, missing-genotype
behavior, and uncertainty after final pedigree consensus. The published
pair-template accuracy remains distinct evidence.

The generator and detector assumptions follow the primary implementation
documentation for [msprime](https://tskit.dev/msprime/docs/stable/ancestry.html),
[Ped-sim](https://github.com/williamslab/ped-sim) and
[hap-IBD](https://github.com/browning-lab/hap-ibd).
