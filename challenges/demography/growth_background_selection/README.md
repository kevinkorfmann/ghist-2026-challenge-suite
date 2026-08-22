# Growth with background selection

Infer growth timing and magnitude while separating directly selected coding sites,
linked background-selection effects, and neutral regions using the coding BED and
recombination map.

The testing VCF contains 42 fully phased diploid samples across chromosome arms
2L, 2R, 3L, and 3R. Of 759,057 records, 738,421 are segregating biallelic sites
after filtering multiallelic and sample-fixed records. Coding intervals cover
18,664,670 bp (17.1% of the supplied sequence). The VCF exposes `AA` and `GT` but
no nonsynonymous annotation field, so coding overlap is treated as an explicit
proxy rather than equating every coding SNP with a nonsynonymous mutation.

The initial signal audit shows strong linked-selection structure: singleton
fraction is 33.4% in coding sequence and 35.2% at zero recombination, while it is
about 25–29% across positive-recombination and noncoding strata. Diversity across
nonempty 100 kb windows correlates 0.63 with local recombination. Primary
inference will therefore use a shared-demography, stratum-aware fit and compare it
against high-recombination noncoding spectra, unstratified spectra, and learned
whole-arm summaries. Whole chromosome arms remain isolated across validation
partitions.
