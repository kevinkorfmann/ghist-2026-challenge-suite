# Bottleneck

Infer the time of the population-size decline and the post-decline population
fraction. Planned comparisons: dadi, moments, momentsLD, simulation-reference
trees, a convolutional context model, and a calibrated ensemble.

The first implemented mechanistic reference is a folded, missingness-projected
moments SFS fit with 24 independent likelihood restarts. Its optimal mutation
scaling estimates ancestral effective size, which converts the fitted `2N` time
to generations. This prevents an otherwise easy but serious unit-conversion error.
