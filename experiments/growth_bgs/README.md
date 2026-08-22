# Growth with background selection evidence

`observed_audit_summary.json` is the compact, versioned intake record for the
testing VCF, coding BED, and Comeron recombination map. It includes source hashes,
record classes, ancestral orientation, coding and recombination strata, and
100 kb diversity diagnostics. The complete window and frequency-spectrum audit
remains on Betty under `artifacts/task3_growth_bgs/`.

The generic VCF audit ran in job `7745884`, the joint resource audit in job
`7745895`, and the compact summary in job `7745897`. The final VCF remains sealed,
and nothing has been submitted.

`stratum_fit_summary.json` records the first mechanistic growth fits from jobs
`7745903` through `7745905`. These independently used 84, 60, and 40 projected
haplotypes, 64 restarts per spectrum, and the exact physical exposure of each
stratum. The whole-genome growth-ratio median is 2.2407, whereas noncoding
positive-recombination quartiles decline from 2.7516 to 1.6011 as recombination
increases. Over the same quartiles, inferred ancestral size rises from 6,682 to
11,945 while inferred modern size remains much more stable. This is diagnostic
evidence of linked-selection bias, not a frozen prediction. Chromosome-arm
consistency and simulation recovery remain required gates.

The refined middle-distance spectrum combines the 1–10 kb and 10–100 kb
noncoding bands, retaining 61,965,066 physical bases while excluding coding,
the most tightly linked 1 kb flanks, and the sparse 2.44 Mb extreme-distance
tail. Jobs `7745927`–`7745929` estimate, respectively, 5,157/1.995 at 84,
5,052/2.010 at 60, and 4,899/2.033 at 40 haplotypes. Their medians are 5,052
generations and 2.010-fold growth. The near-symmetry around 5,000 generations
and close agreement around two-fold growth make 5,000/2.0 the leading testing
hypothesis, but it remains provisional until linked-block and generator-matched
selection recovery pass.
