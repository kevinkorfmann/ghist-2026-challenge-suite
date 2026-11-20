# Task 11 whole-family genotype pipeline smoke, November 20

**Task 11 remains unready.** Two new, deliberately simple synthetic cohorts
successfully traversed genotype preparation, hap-IBD, HBD masking, family
clustering, the frozen dog forest, Bonsai and the existing consensus. The
consensus converted correct sibling calls into incorrect relationships:
**16/96 correct after consensus versus 96/96 before consensus at each HBD
threshold.** This is a concrete composition failure on an easy synthetic case,
not an estimate of challenge accuracy.

## Fixed protocol and independence

The implementation is
[`validate_task11_complete_pipeline.py`](../scripts/validate_task11_complete_pipeline.py).
The Slurm array **10850571** completed both elements successfully, in 43 and 40
seconds, using at most 2.6 GB reported RSS per job. All simulations and estimator
runs occurred on Betty. The isolated source snapshot and large output files are
at `/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-complete-20261120`.
No shared environment or existing snapshot was modified.

Seeds **2026112711** and **2026112811**, eight families per cohort, four observed
full siblings per family, and all other choices were fixed before execution.
Each cohort contains two families with unrelated parents, three with first-cousin
parents and three with sibling parents. Each family's founder haplotypes and
meioses are newly generated and disjoint from every other family. Both complete
cohorts are evaluation only; there is no fitting, parameter search, or row split.

The simulator gene-drops haplotypes jointly through the parent pedigree with
Poisson crossovers on the existing dog average map. Each cohort has 197,855
markers across 1,978.006833142 cM. Marker frequencies are independently drawn
uniformly between 0.05 and 0.5; founder alleles are independent conditional on
each frequency. There is no founder LD. Known true phase and error-free calls
match the phased/error-free input contract; zero missingness is a restricted
subset of the challenge's permitted inputs. Actual phase estimation is not a
required blocker. The frozen result JSON's phrase “no phasing validation”
describes scope only.

The strict repaired preprocessor accepted all markers and imputed zero alleles.
Hap-IBD uses the existing production parameters (1 cM seeds, 0.2 cM extension,
2 cM output, 50 markers, 1,000 bp maximum gap). HBD masking uses the existing
2/5/10 cM thresholds and exposure correction; the original dog forest is loaded
without refitting. Family clustering uses the existing two thresholds and eight
seeds, with eight families supplied as in production. The harness supplies the
VCF sample roster, retaining even samples with no IBD; no truth partition enters
the inference. If clustering is unstable, only the predeclared reference
partition is used diagnostically and strict stability is marked failed.

No final challenge inputs, public challenge answers, or challenge submissions
were accessed. Historical challenge-derived predictions were not evaluated.

## Measured results

| Measure | Seed 2026112711 | Seed 2026112811 |
|---|---:|---:|
| Complete operational stages | Pass | Pass |
| Families recovered exactly | 8/8 | 8/8 |
| Family ARI / minimum clustering stability ARI | 1.0 / 1.0 | 1.0 / 1.0 |
| Missing samples / missed true pairs / false cross-family pairs | 0 / 0 / 0 | 0 / 0 / 0 |
| Forest correct, at each HBD threshold | 48/48 | 48/48 |
| Final consensus correct | 7/48 (14.58%) | 9/48 (18.75%) |
| Wrong final labels with reported confidence zero | 34 | 27 |

The 96 distinct pairs belong to 16 families, not 96 independent simulation
replicates. The three HBD thresholds are repeated measurements of the same
pairs. This sibling-only stress test does not supply balanced eight-label
accuracy or confidence calibration.

HBD detection also ran against actual genotypes. Mean true versus detected HBD
fractions were 0.05121 versus 0.05109 and 0.07102 versus 0.07053 for offspring of
first-cousin parents; 0.20787 versus 0.20644 and 0.25845 versus 0.25690 for
offspring of sibling parents. Both unrelated-parent strata had zero simulated
and detected HBD. Truth HBD lengths are approximated on the 0.01 cM marker grid;
detected totals include only segments at least 2 cM. These favorable values
condition on the simplistic independent-marker founder model.

## Failure mechanism and narrow repair candidates

In `consensus_task11_dog_pedigrees.py`, lines 96–99 accept any permitted Bonsai
relationship tuple unconditionally and replace the forest relationship with
it, even if that label has zero averaged forest support. In these cohorts,
**74 of the 80 wrong final calls are direct accepted-topology overrides**.
Three further errors result from constraining the forest to a topology-implied
meiosis count; the other three result from full-sibling component propagation.
This shows where errors enter the composition. The biological/statistical cause
of Bonsai's incorrect topology has not been established by this smoke.

There is also a deterministic confidence bookkeeping defect: full-sibling
component propagation changes `relation` and `source` at lines 197–199 but does
not recompute `confidence`. Three incorrect propagated `half_sibling` labels
therefore retain **1.0 confidence** from their previous, correct `siblings`
fallback. Recomputing or explicitly clearing confidence whenever a label changes
is a narrow metadata correction; it does not repair relationship accuracy.

The reported HBD stability is 100% in both cohorts despite the large error.
Lines 74–84 apply the same topology override separately to every HBD run before
calculating stability. That statistic cannot independently validate the
topology, and stability must not be presented as accuracy. The existing frozen
forest was itself stable and correct on these cases, making the degradation
particularly easy to trace.

An independent audit should check unsupported topology overrides, constrained
probability support, stale confidence after label changes, and whether the
dog-to-human Bonsai adapter preserves the intended likelihood semantics.
No thresholds or production rules were changed from these outcomes. Any revised
decision policy requires new independent complete-family confirmation, including
heterogeneous relationship labels, before promotion.

## Evidence and remaining blockers

[`pair_before_after.csv`](../experiments/complete_pipeline_2026_11_20/task11/pair_before_after.csv)
records every pair's three forest calls, averaged probability, topology tuple,
final label/source, and reported versus label-specific probability.
[`audit_summary.json`](../experiments/complete_pipeline_2026_11_20/task11/audit_summary.json)
and the per-seed `result.json`, `protocol.json`, `truth.json`, and stage logs
preserve denominators, source/model/map/JAR hashes, and all observed outcomes.
Run [`summarize.py`](../experiments/complete_pipeline_2026_11_20/task11/summarize.py)
to regenerate the before/after audit and verify retained output hashes.

Six focused tests pass for inheritance, disjoint founder families, missing and
merged-family accounting, duplicate rejection and unconditional pair scoring.
Ruff passes for the harness and tests. This validation code does not alter the
existing estimator or historical freezes.

No independent dog founder VCF/BCF was found in the inspected available external
assets. Matched dog founder LD/frequencies, permitted missing-genotype stress,
larger heterogeneous inbred pedigrees covering all accepted labels, and
independent complete-pipeline confirmation remain outstanding. The discovered
consensus degradation is an additional concrete blocker. Supplied family count,
simple sibships and independent-marker founders limit what successful early
stages establish. **No Task 11 final run or submission is justified by this
experiment.**
