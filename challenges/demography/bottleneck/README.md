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
