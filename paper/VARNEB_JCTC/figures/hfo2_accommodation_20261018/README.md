# Local energy release and its geometric accommodation

This revision connects the existing F/A/J energy experiment to the same
geometric T-pattern chart used for the paths. No new electronic calculation,
calculator parameter, model fit or target label is introduced.

## Conclusion and evidence

In the five zero-strain preserving controls, most fixed-cell F-to-A
rearrangement lies outside the three displayed T patterns. The paired
energy lowering is therefore not explained by a single plotted amplitude.
This is a structural observation supporting the local accommodation
interpretation, not a unique microscopic energy partition or a proof of the
origin of the complete strain-dependent barrier.

- Panels(a,b): original F/A/J candidate-energy differences, normalized by
  four formula units per twelve-atom cell. The previous two-panel figure
  remains archived unchanged in hfo2_mechanism_20261018.
- Panels(c,d): the **same** F-to-A displacement projected onto the three
  geometric T patterns and their complete orthogonal complement. These
  lengths are Euclidean norms over the ordered twelve atoms in the original
  T metric, after removing a common translation. They are not per-atom RMS
  displacements. Squared norms, not norms, satisfy the orthogonal sum rule.
- Both rows intentionally use different vertical ranges between channels.
  There are ten physical points and thirty F/A/J geometries, not twenty
  independent points from the two curvature probe steps. No statistics or
  error bars are inferred from distinct coordinate positions.

For preserving F-to-A, the complement norm is0.303110--0.316437Å versus
0.054953--0.062794Å in the triplet. The complement carries95.8909--97.0725%
of the squared displacement; oxygen carries80.8265--81.2416% of the total,
with Hf carrying the rest. Those percentages are **not energy fractions**.
The T F-to-A complement is0.039080--0.065420Å, and its Qz change is
0.020214--0.035555Å, versus less than0.000335Å in preserving.

Only the original two zero-strain centers are covered here. The +0.5%
initialization control and the other strained bottlenecks are separate.
No stable conditional surface, stationary TS, independent predictive gain
or whole-path barrier correction follows from this reanalysis.

## Reproduce

From the repository root, using the existing Python/matplotlib environment:

```text
python -m scripts.analyze_hfo2_local_accommodation --output NEW_DIRECTORY
python -m pytest -q tests/test_hfo2_local_accommodation.py tests/test_hfo2_mechanism.py
```

The output directory must be fresh. `analysis.json` records the source
hashes and all numeric observations; `states.csv` has30rows and `changes.csv`
has20paired contrasts. No native SCF is re-run or re-audited. The original
hash-bound native and optimizer audits are reused, and the reconstructed
F/A/J energies, retained local Q and exact fixed seed cell are checked.

The exports are183×175mm PDF/SVG and300dpi PNG, generated entirely in
Python. Subpanel labels are normal-weight13pt(a)--(d), aligned with the
plot edges; axes align by row and column, top/right spines are visible,
there are no panel titles, and row legends have translucent white framed
backgrounds outside the data. The final PNG was visually inspected: no
label/legend/bar overlap was observed. The original figures and data are
not deleted or rescaled to create this contrast.
