# Sweep readiness follow-up — November 20, 2026

This follow-up improves Tasks 5–8 operational correctness and tests revised
methods on simulations. It does not certify perfect inference, replace current
input provenance, or authorize final-data runs or submission. The earlier
[scientific review](scientific_sweeps_2026_11_20.md) and its failures remain part
of the evidence.

## Provenance and specification limits

The official [sweep specification](https://www.synapse.org/Synapse:syn74803798/wiki/640755)
was retrieved from the public Synapse wiki API, modified October 22, 2026.
Tasks 5 and 6 specify one hard or soft sweep with other loci neutral. Task 7
adds exonic background selection. Only Task 8 explicitly adds bottleneck and
population growth, multiple sweeps, and selected positions restricted to exons.
Testing uses chromosome 7 and mutation rate 1.44e-8; final uses chromosome 9 and
1.29e-8. Twenty-five diploid samples were verified from archived Task 5/6 VCF
headers only; current-file provenance still requires a fresh authenticated run.

During specification review, a grep of cached, single-line wiki JSON
inadvertently displayed embedded public answer sections. No answer values are
recorded here or used for simulation or model choices. Later specification
parsing removed answer sections before display. This report claims a
prespecified protocol, **not investigators unexposed to public answers**.
Targets are drawn mechanically from a separate fixed RNG over accessible
bases, without answer-derived anchors. New confirmation outcomes cannot be
used to tune or select methods.

## Operational changes

The new `python -m scripts.predict_single_sweep` entry point consumes only a
VCF, an accessibility BED and a dataset identity. It runs the full chromosome
scan, chooses the fine-localization region from the scan, writes a provisional
BED and invokes the production submission validator. It has no truth argument.
Existing output artifacts cannot be overwritten accidentally. The candidate
reports input hashes, method identity and the historical estimator-consensus
flag; that flag is an abstention diagnostic, not a calibrated probability.

The chromosome normalization now excludes windows that fail the callable
exposure criterion. Previously these windows could alter every other window's
score even though their own final scores were discarded. Fine localization now
uses the same callable variants as the scan, rejects gaps spanning missing
sequence, and requires sufficient callable exposure on both LD flanks. The old
Task 5 freeze's interval-membership tests now follow BED's `start <= point < end`
convention. Task 7's recovery guard still requires BGS simulations of its exact
estimator, but no longer incorrectly describes bottleneck/growth as a Task 7
specification requirement; the tested demographic history must be explicit.

Example for an explicitly authorized testing run:

```bash
python -m scripts.predict_single_sweep \
  --vcf INPUT.vcf.gz --accessible testing.accessible.bed --dataset testing \
  --bed candidate.bed --out candidate.json
```

This creates a provisional candidate, not a freeze or an upload.

## Prespecified Task 5 confirmation

The [frozen design](../experiments/readiness_2026_11_20/sweeps/task5_design_v4.json)
contains 48 whole-chromosome simulations: two generic demographic histories
(constant diploid N=10,000 and catalog Zigzag_1S14), three selection/age cells,
and eight independent replicates per cell. These are demographic stress
histories, not asserted matches to an undocumented generating history. It uses
25 diploids, the actual HapMapII_GRCh38 map, the official accessibility interval,
mutation rate 1.44e-8 and exact SLiM scale Q=1. Selected positions and simulation
seeds are fixed independently. Simulations condition on selected-allele
frequency at least .95 at sampling; they do not estimate establishment
probability or power for unsuccessful sweeps.

The predictor is locked before biological outputs: 100 kb windows at 20 kb
steps; composite-selected ±750 kb localization; median of composite, largest
callable variant gap and 100 kb flank LD; ±60 kb intervals; estimator span at
most 50 kb for the conditional-consensus diagnostic. Every completed case
passes through the real VCF reader, prediction process, BED writer and output
validator, followed by the independent RAiSD reference. No simulation truth is
passed to either predictor. Scoring occurs only after prediction files exist.

Two initial operational attempts failed before producing biological output:
stdpopsim 0.3.0 could not parse a single-row BED path, and constant-demography
burn-in .1 started too recently for a 3,000-generation mutation. Version 2
passes validated Nx2 intervals; version 3 sets burn-in1. Both failed attempts
and superseded design hashes are retained. Scientific parameters, seeds and
positions were not changed in response to accuracy results. Cross-review then
found that tskit also applies the site-position transform to the declared contig
length. Version 4 preserves exact length and exports binary A/C states to both
predictors, retaining the same genotype information. It asserts equal usable
SNP counts in RAiSD and the production reader. The completed resource-pilot tree
is re-exported rather than resimulated; its accuracy was not inspected before
this correction. The reader now rejects missing/non-diploid calls and unphased
heterozygotes, since dropping them while retaining full physical exposure could
create spurious sweep gaps. Unphased homozygotes remain usable.

Independent cross-review also found that the reader accepted duplicate loci and
symbolic or indel alleles as SNPs. The current reader rejects these ambiguous
inputs. The frozen v4 simulations already export unique A/C SNPs, so this guard
does not change their predictions. A second guard now rereads the raw RAiSD
report, checks finite in-contig peaks, verifies identical SNP counts, and hashes
both report and information files. This prevents stored reference metrics from
silently disagreeing with their raw source. The original v4 runtime is retained;
its reference hashes are attested by the current summary after rereading files,
while future drivers write them at generation time.

A resource pilot precedes wider execution. Any reduction is based on runtime
and memory and retains complete six-cell blocks; missing and failed cases are
reported rather than removed from a favorable denominator. Coverage is reported
both unconditionally and conditional on consensus, alongside abstention,
parameter strata and RAiSD errors. Binomial intervals within exchangeable
strata do not establish narrow error rates with small samples; pooled
heterogeneous strata also receive a conservative independent bounded-variable
uncertainty bound.

The constant-population pilot and next two cases completed in 432.4, 425.8 and
406.2 seconds. The first three Zigzag cases failed during SLiM after 843.1,
960.4 and 723.0 seconds, before producing biological outputs. A live cgroup
trace directly confirms out-of-memory termination for the last case: memory
reached its 11 GiB limit and `oom_kill` increased from zero to one. The other two
terminated with signal 9; their precise cause is not independently established.
A same-seed, unchanged-model retry uses 22 GiB and four CPUs with a two-hour
limit. This pass was then locked to the first six balanced cases, before any
accuracy was inspected; the remaining 42 will not run in this pass. The
[resource record](../experiments/readiness_2026_11_20/sweeps/resource_policy.json)
preserves these decisions before inspecting localization accuracy.

## Completed Task 5 evidence

The [initial-attempt summary](../experiments/readiness_2026_11_20/sweeps/task5_initial_attempt_summary.json)
contains three completed constant-history cases, three failed Zigzag attempts
and 42 cases not run after the resource decision. The completed-case errors are
31.12, 39.89 and 108.04 kb; two of three fixed ±60 kb intervals cover the simulated
target. All three fail the locked 50 kb estimator-consensus criterion, giving
100% abstention and no estimable conditional coverage. All three candidate BEDs
pass the production format validator. The distinct RAiSD reference has errors
142.39, 416.27 and 13.74 kb; one of three is within 60 kb.

There is only one completed replicate per parameter cell, and the pooled
heterogeneous-case Hoeffding interval spans [0, 1]. The original planned
denominator remains 48: observed format-success yield is at least 3/48 and
coverage-success yield at least 2/48. These are incomplete lower bounds, not
claims that unrun cases biologically failed. Even within the six-case resource
block, only half of first attempts produced predictions. Higher-memory retries
are retained separately and do not erase these failed attempts.

The [post-review compatibility run](../experiments/readiness_2026_11_20/sweeps/task5_postreview_verification.json)
passed the same three real VCFs through the current stricter production parser
and reproduced all prediction fields exactly. This establishes implementation
compatibility on those observations; it does not add independent replicates.
The fresh evidence supports retaining Task 5 as provisional, not promoting it to
validated release status.

## Task 6 method development

All 48 earlier standing-variation simulations are explicitly repurposed as
**development data**, including cases once called held out. Their previous
broad recovery failure is not erased. A fixed development grid compares H12
in windows of 50, 100, 200 and 400 SNPs, using all variants or minor-allele
frequency at least .05. This follows the rationale for SNP-delimited
haplotype windows in [Garud et al. (2015)](https://doi.org/10.1371/journal.pgen.1005004),
without assuming Drosophila-calibrated settings transfer to humans. Windows
cannot span inaccessible sequence; windows in which all haplotypes are unique
are marked uninformative. Tied maxima and their deterministic selection rule
are recorded.

Candidate selection uses median localization error plus one quarter of the
90th-percentile error on these development simulations. Any selected candidate
still requires new exact-scale confirmation. Q10 development results cannot
certify Q1 accuracy or chromosome-scale transfer. Tasks 7 and 8 retain the
previously measured BGS/multiple-sweep limitations; this follow-up does not
convert small-segment results into full-chromosome calibration.

## Completed Task 6 development evidence

All 48 old simulations completed all eight new SNP-window comparisons. The
best development objective selected 50 SNPs without the MAF filter, but its
median error was 165.18 kb, 90th percentile 1.853 Mb, and only 19/48 intervals
within ±75 kb would cover the selected site. It is not promoted. The existing
median-three estimator on the same development simulations had median 37.05 kb,
90th percentile 2.437 Mb and 30/48 within 75 kb. The new candidate trades a
smaller extreme tail for substantially worse typical precision; this is not
credible evidence of an overall operational improvement.

Looser end-frequency conditioning is associated with failure: among twelve
simulations with minimum end-frequency thresholds .05, .2 or .4, only two
median-three estimates were
within 75 kb. However, even the 16 simulations conditioned at .95 had a 1.430 Mb
90th-percentile error. These are conditioning thresholds, not reconstructed allele trajectories.
These posthoc, correlated parameter strata do not prove a causal effect. Only 4/48 passed the pre-existing 50 kb consensus diagnostic;
all four were within 75 kb, with 91.7% abstention. Four selected successes cannot
certify reliable conditional coverage.

The Task 6 observed-transfer script also contained a real scale mismatch:
it compared observed 25 kb-window scores with simulated 100 kb-window scores.
It now requires explicit matching scale metadata, computes ensemble errors
correctly, treats empty matching strata as failed gates, and does not describe
retrospective/development records as untouched validation. Historical artifacts
are preserved; this repair does not turn earlier failed recovery into a pass.

Cross-review found that the first development summary trusted cached errors and
omitted the coordinate origin needed to recompute them. Revised records retain
the segment origin, length, source-record hash and tree hash. Aggregation
recomputes every error from the relative prediction and absolute truth, rejects
nonfinite or negative values, and requires identical candidate settings across
replicates. The same 48 development trees were rerun into separate
[v2 evidence](../experiments/readiness_2026_11_20/sweeps/task6_development_v2.json).
All 384 candidate predictions and errors exactly match the archived original
records; the scientific conclusion is unchanged.

## Reproduction and code checks

The [verification record](../experiments/readiness_2026_11_20/sweeps/verification.json)
lists all attempted scheduler jobs, including setup failures. The exact
[source snapshot](../experiments/readiness_2026_11_20/sweeps/source_snapshot_v4.tar.gz)
and [software hashes](../experiments/readiness_2026_11_20/sweeps/software.json)
preserve the confirmation runtime independently of concurrent edits to shared
validators. Extract the snapshot in an isolated directory, set its `src` and
root on `PYTHONPATH`, and use the committed Slurm templates with the matching
frozen design. Heavy simulation and reference scans run through Slurm only.

Forty-one focused sweep tests pass, including genotype/coordinate round trips,
callable-mask invariance, source/identity guard behavior, half-open interval
semantics, matching Task 6 score windows, empty-stratum handling, duplicate and
non-SNP input rejection, and corrupt reference/development-metric rejection. Sweep Ruff
checks and Slurm template syntax checks pass. These checks establish
software behavior; biological recovery is reported separately.
