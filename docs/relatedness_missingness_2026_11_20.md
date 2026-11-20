# Relatedness missing-genotype diagnostic — November 20, 2026

The official relatedness input contract allows missing genotypes in otherwise
phased, error-free VCFs. The shared preprocessing repair rejects an implicit
imputation decision; this experiment measures two explicit policies on the
already generated synthetic Task 11 cohorts. It is additional sensitivity
analysis of those same families, not new biological confirmation.

Before observing these outcomes, the protocol fixes whole-genotype dropout
probabilities of 0.005 and 0.02, independently across markers and individuals.
Both original cohorts (seeds 2026112711 and 2026112811) are included. The masking
seed is the cohort seed plus 1,000,000 plus the index of the dropout rate. Each
masked VCF is analyzed under both `drop-site` (remove incomplete markers) and
`major` (the explicit historical major-allele replacement). No method or numeric
threshold is selected from these results.

The two branches use the same masked calls, frozen dog model and map, hap-IBD
parameters, eight-family clustering setup, and 5 cM HBD masking. Source, input,
baseline result and artifact hashes are recorded before execution. Any failure
is retained and missing true pairs remain in the accuracy denominator. The
already failing pedigree consensus is not reassessed here; this diagnostic
isolates genotype preparation, detection, clustering and pair classification.

All 32 sampled individuals per cohort are retained. Dropping a marker whenever
any sample is missing has expected retained fraction `(1-p)^32` under this
masking process. Retention will be much lower for larger cohorts at the same
per-individual missingness rate. Missingness correlated with genotype, position,
family or individual is not covered. Neither these founder markers without LD
nor the sibling-only pedigrees establish complete Task 11 accuracy or transfer
to Tasks 10 and 12.

The reproducible harness is `scripts/validate_relatedness_missingness.py`;
the Slurm array is recorded in
`experiments/complete_pipeline_2026_11_20/missingness/run.sbatch`. Its four array
elements compare two policies each, reusing two independent original cohorts.
The local masking/preparation tests pass (18 combined tests). Execution uses a
separate remote source snapshot and preserves the original whole-family runs.

Slurm array 10850656 completed all four jobs in 48–53 seconds each (exit 0).
All eight policy executions recovered all eight families and classified all
48 sibling pairs correctly in each original cohort. No pair was missing and
no cross-family pair was introduced. These repeat the same 96 distinct sibling
pairs; they are not 384 or 768 independent biological observations.

| Missing genotypes | Drop-site marker retention, two cohorts | Major-replacement marker retention | Pair accuracy, either policy |
|---|---|---|---|
| 0.5% probability | 85.14%, 85.11% | 100%, 100% | 48/48 in each cohort |
| 2% probability | 52.11%, 52.54% | 100%, 100% | 48/48 in each cohort |

The [compact result](../experiments/complete_pipeline_2026_11_20/missingness/summary.json)
and per-run records preserve the masks, settings and failed-stage policy. This
narrow sibling diagnostic does not distinguish the policies by relationship
accuracy and does not support switching the production default. The large loss
of markers under drop-site is real, even when this particular outcome stays
correct. Higher missingness, larger cohorts and distant relationships still
need testing. No final challenge data or competition submission is involved,
and no pipeline is promoted by this test.
