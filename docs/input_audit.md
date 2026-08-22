# Input and evaluation audit

Before fitting any model, every dataset receives a frozen audit. Findings become
simulation factors or validation strata and are not used as answer-specific cues.

## VCF checks

- sample count, names, ordering, ploidy, and phased-call fraction;
- partial versus complete missingness by sample and genomic position;
- biallelic SNP, multiallelic, indel, and complex-record fractions;
- FILTER, INFO, FORMAT, annotation, and ancestral-state fields;
- contig names, header lengths, covered coordinate range, and discontinuities;
- intervariant-spacing and folded-frequency distributions;
- duplicated positions, non-monotone records, and coordinate normalization;
- invariant/inaccessible intervals versus biological low-diversity intervals;
- coding and recombination-map coordinate conventions and assembly compatibility.

## Generator and split checks

- fixed versus randomized ancestral size and nuisance parameters;
- mutation and recombination-rate misspecification;
- simulator, ancestry model, gene-conversion, and selection-model mismatch;
- chromosome-specific accessibility and recombination distribution shifts;
- sample-label or file-metadata cues that disappear on the final dataset;
- independent seeds and parameter cells across training and validation;
- pedigree/family separation for relatedness tasks;
- one-based VCF positions versus zero-based half-open BED intervals.

## Evaluation checks

- exact parameter keys, units, rounding, and schema;
- point-locus answers versus interval-overlap scoring for sweeps;
- handling of ties, duplicate pair labels, and relationship symmetry;
- sensitivity to interval padding and threshold selection on held-out simulations;
- checksum and immutable provenance for every generated prediction.

