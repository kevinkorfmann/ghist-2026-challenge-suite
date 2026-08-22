# Bottleneck

Infer the time of the population-size decline and the post-decline population
fraction. Planned comparisons: dadi, moments, momentsLD, simulation-reference
trees, a convolutional context model, and a calibrated ensemble.

The first implemented mechanistic reference is a folded, missingness-projected
moments SFS fit with 24 independent likelihood restarts. Its optimal mutation
scaling estimates ancestral effective size, which converts the fitted `2N` time
to generations. This prevents an otherwise easy but serious unit-conversion error.

Linkage summaries use distance-stratified streaming reservoirs over true local
variant pairs from 100 bp to 1 Mb. The identical sampler is applied to observed
and simulated records, avoiding short-range linkage loss from genome-wide SNP
subsampling.

Half of the training simulations remain exactly on the announced mutation and
recombination rates. The other half form labeled robustness strata spanning
global rate error, piecewise recombination heterogeneity, missing calls, and
smooth low-diversity tracts. Validation reports the nominal and perturbed strata
separately so robustness cannot conceal a loss of nominal precision.

Observed inference preserves the simulation scale: the 100 Mb chromosome is
represented by fifty contiguous 2 Mb blocks, each with 64 local windows. Every
nonempty block contributes to a robust precision aggregate and is retained in a
diagnostic artifact. The whole chromosome is also fit directly by the projected
SFS reference, leaving both local and global information available.

The testing VCF contains 57 diploid individuals, 1,194,350 fully phased and
complete biallelic SNP records, and correct ancestral-allele annotations. The
downloaded header has one malformed metadata line and 7,715 adjacent repeated
coordinates. A scheduled normalization repairs only that header and proves the
SHA-256 of every record byte is unchanged. Summary and FastCXT inputs retain the
original coordinates. Genealogy inputs retain every record while assigning
stable sub-base order inside repeated-coordinate groups.
