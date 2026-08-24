# Single hard-sweep evidence

The official testing data represent one completed hard sweep on GRCh38
chromosome 7 in a population with a bottleneck and growth. The VCF contains 25
phased diploids and 702,459 records across 159,345,973 bp. Of these, 701,382
fully called biallelic sites are retained and 1,077 multiallelic records are
excluded. The chromosome 7 accessibility mask begins at 35,326 bp and extends
to the chromosome end. The final chromosome 9 VCF remains sealed, and nothing
has been submitted.

The multiscale scan combines five distinct signals rather than allowing one
statistic to determine the answer: nucleotide diversity, segregating-site
density, folded Tajima's D, phased haplotype homozygosity, and pairwise linkage
disequilibrium. Windows from 25 kb through 500 kb converge on one region near
38.2–38.4 Mb. At 50 kb, the leading window is 38.200–38.250 Mb with diversity
of 0.00001997, Tajima's D of -2.65, and H12 of 0.709. At 25 kb, H12 reaches
0.961 in the same local trough. The next chromosome-wide broad-scale peaks are
substantially weaker.

A split-flank LD scan and largest-variant-gap calculation refine the region.
The 100 kb flank statistic peaks at 38,224,000 bp, while the largest local gap
spans 38,245,545–38,268,357 bp and has midpoint 38,256,951 bp. Wider LD flanks
shift left and are retained as evidence of heterogeneous recombination rather
than averaged away.

The independent RAiSD 2.9 reference scan evaluates 701,333 locations with its
default 50-SNP windows. Its maximum is 38,240,723 bp with a composite score of
204.9. The next 100 kb-separated peak has score 62.74, and the strongest peak
outside the primary region has score 10.26. The reference maximum is 16 kb
from the largest-gap midpoint and 17 kb from the fine split-flank maximum.

The neutral chromosome background, excluding 35–41 Mb, has diversity
0.0008163, singleton fraction 0.2947, and intermediate-frequency fraction
0.2331. Among nine published human histories simulated under the species-wide
mutation and recombination rates, `Zigzag_1S14` is the closest by folded-spectrum
and diversity distance. This background anchors selection recovery while
population history, recombination, sweep age, strength, and completion are
varied as nuisance parameters. A chromosome-7 recovery campaign using the
HapMapII GRCh38 map selected the median of the gap, composite, and split-flank
estimates without access to its ten evaluation cells. On those untouched cells
it has 7.1 kb median absolute error, 45.2 kb 90th-percentile error, and 90%
recovery within 50 kb. An exploratory scale-factor campaign exposed long-tail
errors in weak synthetic sweeps instead of passing its deliberately strict gate.
A 192-replicate, unaccelerated chromosome-map campaign then separated 64
calibration, 64 pilot, and 64 untouched evaluation cells. On the untouched
tranche, overall median error is 27.0 kb and the 90th percentile is 114.8 kb.
Twelve cells match the observed estimator-agreement pattern; all twelve recover
within 60 kb, for a mean challenge score of 0.625. Every predeclared gate passes.

Interval width is optimized against the published GHIST sweep evaluator rather
than chosen visually. Its length penalty selects a 25 kb half-width on the
calibration cells. A predeclared consensus guard first expands that to 35 kb so
the candidate contains the observed localizers. Independent exact-scale
calibration selects a 60 kb half-width. The frozen testing interval is therefore
chr7:38,196,951–38,316,951, centered at 38,256,951 bp.

The validation strategy follows the multi-signature logic of
[RAiSD](https://doi.org/10.1038/s42003-018-0085-8), the spatial-statistic design
of [diploS/HIC](https://doi.org/10.1093/molbev/msy006), and the nested-window,
confounder-aware principles of
[Flex-sweep](https://doi.org/10.1093/molbev/msad139). The forward-time matched
selection smoke test, parameter-cell-isolated localization campaign, and
untouched exact-scale gate have passed. The method and BED interval are frozen
in `task5_freeze.json`; nothing has been submitted.
