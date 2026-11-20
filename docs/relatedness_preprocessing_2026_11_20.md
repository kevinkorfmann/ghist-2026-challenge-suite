# Relatedness genotype preparation repair — November 20, 2026

The shared Tasks 10–12 preparation script previously replaced `/` with `|` for
every genotype. An unphased heterozygote such as `0/1` therefore became an asserted
haplotype phase without inference. It also replaced missing alleles with the
site's major allele without an explicit request. Both behaviors can change the
haplotype evidence supplied to the relationship estimator. This finding is a
code-path defect; it does not establish that a particular archived input
contained unphased heterozygotes or that every historical output was affected.

The repaired `scripts/prepare_hapibd_vcf.py` rejects unphased heterozygotes and,
by default, all missing alleles. Homozygous calls may safely be normalized from
`0/0` to `0|0` because both haplotypes carry the same allele. An explicit
`--missing-policy major` option retains the historical major-allele replacement
for reproducibility, reports replacement and skipped-site counts, and still
refuses to invent phase for heterozygotes. That option is not a validated
imputation method. Sites without any observed alleles are skipped only in this
explicit mode or the separate `--missing-policy drop-site` mode, which removes
an entire marker if any sample has missing data and never imputes alleles. A
file with no usable records is rejected. Dropping markers can sharply reduce
coverage in large cohorts; it also requires independent accuracy assessment.

The [current relatedness specification](https://www.synapse.org/Synapse:syn74803798/wiki/640950),
checked through its public API, describes phased, error-free VCFs while allowing
missing genotypes. Thus supplied true phase matches the challenge contract;
phase estimation is not a required challenge stage. Missingness handling is a
real remaining validation requirement. The strict default deliberately refuses
to select an unvalidated missing-data policy implicitly.

The script validates diploid genotype fields and sample headers, preserves
non-genotype FORMAT values, supports plain and gzip-compressed VCFs, and only
replaces an output after the entire input passes. A malformed later record
cannot leave a partial result or overwrite the previous output. The human and
dog Slurm wrappers use the default strict policy and save preparation diagnostics.
Their output is named `input.phased.vcf.gz`; historical files are preserved.

The behavior follows the [official hap-IBD input requirements](https://github.com/browning-lab/hap-ibd#required-parameters):
input must have phased genotypes without missing alleles. Checking delimiters
cannot establish that a caller's reported phase is biologically correct. This
repair supplies no phasing algorithm and no new relationship-accuracy claim.

Validation: 15 focused tests passed locally, including phase preservation,
rejection under both missing-data policies, explicit imputation, all-missing
sites, malformed genotype/header contracts, plain/gzip input, and atomic output
preservation. Ruff and shell syntax checks passed. The new Task 10 and Task 11
synthetic whole-family workflows use this repaired preparation path.

No current final data, final inference, upload, or competition submission was
used for this repair.

Independent cross-review subsequently found that explicit independent `PS`
blocks were accepted as chromosome-wide haplotypes. The
[VCF specification](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.tex)
defines phase orientation within a phase set, while
[hap-IBD's GT parser](https://github.com/browning-lab/hap-ibd/blob/master/src/vcf/VcfRecGTParser.java)
does not connect separate `PS` blocks. The current preparation source therefore
requires a single explicit or implicit phase set per sample and chromosome
across retained heterozygous calls. Mixed implicit/explicit sets and multiple
explicit sets fail before publishing output. Homozygous sites do not constrain
orientation; dropped sites do not connect retained blocks. Ordinary phased VCFs
without `PS`, one consistent `PS` per sample/chromosome, and differing sets on
different chromosomes remain supported. This validates an input contract and
does not estimate or certify biological phase.

The follow-up combined suite passed **170 tests**. Tests include both block
orders, absent/missing `PS`, independent chromosomes and samples, homozygotes,
missing-site removal, explicit major-allele imputation, invalid informative
phase-set IDs and preservation of previous output on failure. The original
experiments and their runtime source snapshots remain unchanged; they did not
run this later phase-set guard. Before/after hashes and preserved-file checks
are recorded in the
[cross-review verification](../experiments/complete_pipeline_2026_11_20/cross_review/verification.json).
