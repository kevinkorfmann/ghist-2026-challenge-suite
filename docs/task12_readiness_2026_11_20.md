# Task 12 independent-reference repair, 20 November 2026

The corrected relationship-model candidate failed its prespecified improvement
gate and is not promoted. Exact-label accuracy was 37.37%, compared with 37.64%
for the historical model under the same independent correction reference. This
fixed comparison does not validate the full genotype-to-family pipeline.

## Prespecified experiment

The [protocol](../experiments/readiness_2026_11_20/task12/protocol.json) was frozen
before execution. It locks the two original model artifacts, six original
calibration segment files, source code, pedigree definition, Ped-sim binary,
sex-specific map and interference parameters by SHA-256. Training uses only
calibration replicates 0–5. Each augmented training replicate is corrected using
the other five replicates; an evaluation replicate is never its own clean
reference. Row-random early stopping is disabled. The historical direct
histogram-gradient classifier architecture, 240 iterations, class balancing,
feature set and blend weight 1.0 are fixed without a model or weight search.

The frozen empirical background-bin profile is inherited from the old artifact.
It was originally estimated from unlabeled challenge and reference cohorts; this
repair reads no new challenge inputs and uses no public answers. This inherited
profile is an empirical approximation, not a population-genetic model of shared
founder ancestry.

Only after fitting, Ped-sim generates new interference-model replicates with
seeds 2026114011–2026114014. Each seed is tested under the five previously defined
final contamination regimes. All regimes for a seed share its pedigree segments;
the twenty comparisons are not twenty independent experiments. The independent
correction reference is original calibration replicate 0 for every holdout.
Both old and new models receive identical corrected inputs, and both use the
production probability-blending/prediction helper. Every resulting CSV is
validated against the complete pair set in that simulated segment input.

The prespecified improvement gate requires a mean exact-label gain of at least
2 percentage points, positive mean gain in at least 80% of independent seeds,
no regime losing more than 1 percentage point, and no mean meiosis-accuracy loss
over 0.5 percentage points. Passing this relative gate alone would not establish
high absolute accuracy or justify release. Exact and balanced accuracy, meiosis
accuracy, class recall, Brier score, negative log likelihood and confidence
calibration error are reported without fitting a calibrator to these holdouts.

## Execution and results

Job 10850487 completed with exit code 0 in 2 minutes 17 seconds on eight CPUs,
using 4,075,964 KiB peak RSS (3.89 GiB). Its scheduler-only time limit was shortened
from 30 to 15 minutes to permit backfill; scientific settings and source hashes
were unchanged. Training took 54.11 seconds for 203,950 training rows and retained
the historical estimator random state 122228; augmentation uses seed 122226.

The four new replicates contained 27,212 distinct input pairs. All 40 output CSVs
(two models × four seeds × five regimes) passed schema and supplied-input pair
coverage validation. Mean results, weighting each seed and regime equally, were:

| Metric | Historical model | Corrected refit |
|---|---:|---:|
| Exact-label accuracy | 37.636% | 37.372% |
| Balanced exact-label accuracy | 37.694% | 37.372% |
| Meiosis accuracy | 53.638% | 53.533% |
| Confidence ECE, 10 equal-width bins | 0.0241 | 0.0328 |
| Multiclass Brier score, summed over labels | 0.6697 | 0.6720 |
| Negative log likelihood | 1.4968 | 1.5046 |

The mean exact-label change was −0.264 percentage points. All four independent
seeds had negative mean changes (−0.179, −0.209, −0.409 and −0.258 percentage
points); the worst individual regime lost 1.189 percentage points. Three of four
improvement criteria failed. The independent-reference repair corrects the
experimental design, but this fixed refit does not solve the accuracy problem.
These results do not support calling either model excellent or release-ready.
Confidence statistics are descriptive and do not certify calibrated probabilities.

The [complete comparison](../experiments/readiness_2026_11_20/task12/holdout.json),
[aggregate](../experiments/readiness_2026_11_20/task12/aggregate.json) and
[training provenance](../experiments/readiness_2026_11_20/task12/training.json)
are preserved. Run `python3 experiments/readiness_2026_11_20/task12/aggregate.py`
to reproduce the aggregate from the recorded holdout report. None of these
holdout results was used to change the fitted model or settings.

The snapshot is
`/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-task12-20261120`.
Its candidate is `artifacts/task12_reference_v2/model.joblib`, SHA-256
`84c573b49a0230d5937f5557d0f2780c1f5f0dbafeb614b7dd55360a05a0b1dc`.
No production selection, historical freeze, challenge prediction or competition
submission is changed.

The [job script](../experiments/readiness_2026_11_20/task12/run.sbatch) is the exact
reproduction command sequence. From its prepared snapshot, the original launch was:

```sh
source /etc/profile
sbatch experiments/readiness_2026_11_20/task12/run.sbatch
```

A different output snapshot is required for a repeat; preserve the original
candidate, frozen protocol and report. A path-adjusted job script also requires
a new protocol copy recording that script's updated hash; numerical settings and
seeds must remain fixed for reproduction. Repeated seeds are reproductions, not
new independent confirmation. The historical artifacts and external
simulator inputs are resolved under the explicit `TASK_ROOT` in the job script.
The interpreter is the original Betty environment, whose package versions are
also recorded for the [Task 9 runtime](../requirements/task9-release-2026-11-20.txt).

## Remaining full-pipeline work

The available independent relatedness assets contain exact pedigree IBD segments
and genetic maps. They do not contain independent founder genotypes or a complete
joint family/population simulation. Existing relatedness VCF artifacts are
challenge-derived and cannot serve as untouched synthetic confirmation.

The missing benchmark must generate correlated founder background and a complete
pedigree, with truth for all required cohort pairs. It must then run the actual
VCF preparation, hap-IBD detection, segment/shape summaries, background correction,
relationship model, branch/Bonsai reconstruction and degree-guarded consensus,
followed by the independent VCF-sample pair-coverage validator. Phasing and marker
ascertainment must match the assumed input contract or be explicit perturbations.
Simply adding independent extra segments to each relationship pair cannot test
these dependencies or family reconstruction. This benchmark requires new data
and a prespecified simulation design; it is outside this bounded component refit.
