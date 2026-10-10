# Evaluate frozen forecasts without manufacturing an advantage

`vcneb.prediction_evaluation` implements ordinary deterministic error
bookkeeping for the B0/B1 and strong-reference B2–B5 comparisons registered
in [protocol v1](HFO2_PREDICTION_PROTOCOL.md) and [v2](HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md).
It is calculator/backend independent. It neither supplies the unfinished
material models nor freezes them, reads a holdout, or constitutes new theory.

## Inputs are a complete panel, not only successful predictions

Each `FrozenChannelForecast` identifies a registered model and case,
one unchanged model-artifact SHA256, the shared comparison-contract hash,
and a numerical prediction or an explicit abstention. Preserve a negative
signed prediction with its refusal reason, not a clipped zero barrier.
Different Cmma registrations are different model IDs; do not silently
merge their subspaces or select a favourable one from unseen labels.

Each `AuditedBarrierReference` supplies a directional barrier in eV/f.u.,
a measured deterministic error bound in the same units, the same comparison
contract and the source-audit hash. **These declarations do not audit DFT.**
The upstream study adapter must independently check the physical inputs,
same mechanical family/condition, initial state, formula normalization,
phase, sampling and numerical error. A force threshold is not a barrier
error. For finite pressure, use audited H=E+PV consistently.
The numeric constructor cannot detect that a caller copied a force value
into an energy-bound field; this provenance/unit misuse must be rejected
by the upstream material audit, not inferred to be checked by this scorer.

Every registered case needs every model's record, including failures.
Duplicate/extra/missing cases or models, mixed comparison hashes and changed
model artifacts are rejected. Predictions outside a model's valid domain
must already carry the frozen abstention reason. The scorer does not
retrain, repair or reinterpret such a failure after looking at a reference.

## One shared unknown DFT reference

For a deterministic reference interval y in[L,U] and fixed prediction p,

```text
absolute-error range = [distance(p,[L,U]), max(|p-L|, |p-U|)]
gain(candidate,baseline;y) = |baseline-y| - |candidate-y|
```

The gain uses the **same y** in both errors. It is incorrect to assign
independent reference fluctuations to the two comparisons. With identical
predictions, the shared-reference gain is exactly zero even when each error
separately has a nonzero interval. The gain extrema occur at endpoints and
piecewise-linear breakpoints; this implementation evaluates them explicitly.

For several valid registered baselines, also score

```text
min_baseline |baseline-y| - |candidate-y|
```

The nearest-baseline envelope changes at adjacent prediction midpoints,
which are included as breakpoints. This conservative **scoring oracle** is
not a post-holdout deployable model selector. Keep every pairwise comparison
and model's coverage/error alongside it, rather than reporting only a weak
fixed reference. An unavailable strong baseline is not a poorly performing
baseline: comparisons with it remain unavailable.

## Report coverage and limits

The result gives per-case raw absolute errors and reference-error ranges,
per-model coverage and maximum absolute error, pairwise gains and a
best-available-baseline envelope. Different-coverage error maxima are not a
fair accuracy comparison. Abstentions have null error, not zero error.

The deliberately strong `full_panel_strict_advantage` flag requires all
registered models/cases available and positive gain over the best baseline
throughout every reference interval. Its default 1e-12eV/f.u. arithmetic
tolerance only prevents floating-roundoff claims; it is not DFT uncertainty,
a changed δ=0 H1 criterion or a prospective physical error guarantee.
Passing this numerical flag alone still does **not** establish equal
information/cost, true independent label blinding, mechanism, generalization,
or JCTC readiness. A smaller centre-point error with unresolved gain is
reported but not called a resolved quantitative advantage.

Costs stay in [the separate recorded-cost contract](NESTED_RESPONSE_CONTROLS.md),
including shared, unused and failed preparations. No model gains a free
curvature matrix just because this scoring function itself makes no DFT calls.

## Usage and implementation evidence

```python
from vcneb.prediction_evaluation import (
    AuditedBarrierReference, FrozenChannelForecast, evaluate_frozen_panel,
)
# Construct records only from actual frozen models and audited references.
report = evaluate_frozen_panel(forecasts, references,
    candidate_model="B5_path_adaptive",
    baseline_models=["B0", "B1_fixed", "B1_follow", "B2_T",
                     "B2_Cmma_gauge1", "B2_Cmma_gauge2",
                     "B2_Cmma_gauge3", "B2_Cmma_gauge4", "B3", "B4"])
```

No real held-out HfO2 labels are used by the bounded implementation check:

```sh
python -m scripts.check_prediction_evaluation --output /new/scoring-check.json
python -m pytest tests/test_prediction_evaluation.py -q
```

It checks48random shared-reference intervals against independent10001-point
grids, an additional non-grid-aligned cusp, the envelope midpoint extrema,
zero gain for identical predictions,
strong-reference/accurate-B0 rejection of a manufactured advantage, coverage
failures, complete-panel gates and finite arithmetic. Grid-extremum distances
are discretization checks, **not** DFT error bars or material model accuracy.
The illustrative T/Cmma labels and hashes are synthetic placeholders, not
computed reference models or actual forecast freezes. The existing full-
joint harmonic null must still be respected in actual material analysis.

Material B2–B5 training/selection, the complete frozen batch, subsequent
independent labels and the prospective comparison remain required. This
interface does not authorize any held-out task or change the finite budget.
