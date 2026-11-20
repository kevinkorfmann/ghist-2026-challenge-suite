# Parallel relatedness follow-up — November 20, 2026

Six agents were assigned distinct work after checking the active challenge
review. At most three child agents ran concurrently. The existing review
continued Tasks 1–8 and 12; this follow-up covered complete-family experiments,
the Task 9 handoff, and shared relatedness correctness. This report supplements
[the challenge readiness table](readiness_2026_11_20.md).

| Workstream | Measured outcome | Readiness effect |
|---|---|---|
| Task 9 offline package and independent package review | 64 package/release tests pass; versioned model, map, source, runtime and evidence identities verified | Validated candidate packaged; current final-data prediction remains pending |
| Task 10 complete-family pipeline | One independently generated 18-person family; all 153 pairs produced; 58/153 exact labels and 117/153 meioses correct | Execution demonstrated, relationship accuracy inadequate |
| Task 11 complete-family pipeline | Two independently generated cohorts, 16 families and 96 distinct sibling pairs; all families recovered and all pair-model calls correct; final consensus only 16/96 correct | Consensus fails; sibling-only diagnostic does not certify the broader task |
| Task 11 consensus audit | Coverage/probability checks and confidence reporting repaired; one historical cohort now correctly refuses a zero-support topology constraint | Freeze requires independent complete-task evidence; separate development policy labels only 7/96 pairs and abstains on 89 |
| Shared haplotype-to-IBD conversion | Correct matching reduced overcalled IBD2 in 72/96 pairs; frozen sibling labels stayed unchanged | Corrected feature semantics require new downstream validation; topology/consensus were not rerun |

The coordinator also repaired shared genotype preparation and measured two
explicit missing-call policies. Independent whole-genotype dropout at 0.5% and
2% was applied to the two existing Task 11 cohorts. All eight policy executions
retained family clustering and sibling calls. Dropping incomplete sites retained
about 85% and 52% of markers, respectively; major-allele replacement retained
all sites. These are repeated measurements of the same 96 sibling pairs, not
new independent biological evidence. Neither policy was selected from the
results, and the production default rejects missing calls until an explicit,
validated policy is chosen.

The official genotype contract supplies phased, error-free VCFs but permits
missing genotypes. The preparation repair preserves phased calls, rejects
unphased heterozygotes and reports explicit missing-call handling. It does not
manufacture phase by replacing genotype separators. See the
[preprocessing report](relatedness_preprocessing_2026_11_20.md) and
[missingness experiment](relatedness_missingness_2026_11_20.md).

The Task 10 and Task 11 measurements preserve their original source snapshots.
They do not validate subsequently repaired converter or consensus code. Task
10's low accuracy already occurs in the flat classifier before the pedigree
converter, so the converter bug cannot explain that entire failure. Task 11's
original cohorts become development evidence once used to diagnose or change
the policy; replaying them is not fresh confirmation.

Detailed handoffs and evidence:

- [Task 9 reproducible package and deferred final-run handoff](task9_submission_handoff_2026_11_20.md)
- [Task 10 complete-family result](task10_complete_pipeline_2026_11_20.md)
- [Task 11 complete-family result](task11_complete_pipeline_2026_11_20.md)
- [Task 11 consensus audit and development-policy limits](task11_consensus_audit_2026_11_20.md)
- [Shared IBD conversion repair and measured audit](haplotype_ibd_conversion_2026_11_20.md)

The full workspace regression suite passed **909 tests** (two existing PyTorch
nested-tensor warnings). The final consensus/freeze checks additionally passed
48 focused tests after that run. Ruff passed for all files changed by these
relatedness workstreams; both changed Slurm wrappers passed shell syntax checks.
The Task 9 package was reverified against the trusted manifest digest, and all
four checksums in its transfer manifest passed. These software checks do not
establish additional scientific accuracy.

No current final challenge inputs were downloaded or used for inference, and no
competition submission was made in this follow-up. Shared uncommitted work from
the other active review was preserved.
