# Large family with high background relatedness

Infer pedigree relationships while distinguishing recent pedigree IBD from elevated
population-level background sharing.

The historical method estimates the unlabeled background length spectrum from the
cohort, trains on short-segment contamination, and uses pedigree topology only
when it agrees with the pairwise meiosis degree. Its original reported improvement
used the clean evaluation replicate as the correction reference; that comparison
does not establish performance with an independent reference.

The [scientific audit](../../../docs/scientific_relatedness_2026_11_20.md) documents
this problem. A fixed, separately versioned refit on four new pedigree seeds
failed its improvement gate: exact accuracy was 37.37%, versus 37.64% for the
historical model under the same independent correction reference. See the
[repair and validation report](../../../docs/task12_readiness_2026_11_20.md).
This is relationship-model component validation. The complete VCF, IBD detection,
and family-topology pipeline still needs an independent end-to-end benchmark.
Historical freezes are preserved, and no replacement model is promoted here.
