# E165 — does the original local quadratic describe the released DFT candidates?

Retrospective agreement check on already seen training data. Previous goal turn
E164 was PROGRESS (atomic/cell signatures, e9c1de67). This analysis reuses the
existing conditional-quadratic API and native audit records; it does not fit
coefficients, select a new reference, generate a held-out label or call DFT.

## Result and next decision

The original first-step retained axes are Hessian eigendirections. Thus Q-R
coupling is numerically zero at that original center, and its Schur curvature
correction is approximately 1e-28 eV/Angstrom squared in norm. Releasing R in
this formal quadratic changes the constant energy offset, but cannot supply
a substantial new shape correction in that eigenplane. This is a construction
property, not a newly demonstrated prediction improvement.

After subtracting each model's and observed candidate's seen Q0 energy, the
shape still disagrees with the saved DFT candidates. The table uses all six
common numeric-Q points, not points selected by error. T and preserving axes
have different physical anchors; equal numeric Q is not the same structure.
All nine available T and six preserving candidates remain in analysis.json.

| Center | Curvature probe step, Angstrom | Observed Q0 release | Formal quadratic Q0 release | Maximum absolute zero-aligned shape error |
|---|---:|---:|---:|---:|
| T | 0.01 | -0.652142 | -0.866958 | 8.524106 |
| T | 0.02 | -0.652142 | -0.874429 | 8.587800 |
| Preserving | 0.01 | -49.484930 | -54.521222 | 20.521839 |
| Preserving | 0.02 | -49.484930 | -54.705828 | 20.720640 |

All three energy columns are meV per twelve-atom cell, not per HfO2 formula
unit. The largest common-cohort shape discrepancies occur at Q=(-0.2,0).
Zero alignment uses the already seen relaxed origin and is only a shape
diagnostic, not a new calibration-free forecast. The model's original nonzero
released gradient is retained in both its offset and effective retained gradient.

Observed axial candidate-energy second differences at displacement 0.1 Angstrom
are approximately (-0.405128,-0.243087) eV/Angstrom squared for T and
(-0.337763,-0.184576) for preserving. These are finite candidate differences,
not infinitesimal Hessians at certified conditional minima. Raw retained
partial-gradient secant matrices are also retained, without symmetrizing them
into a certified Hessian. Finite released residuals prevent assuming exact
envelope derivatives.

Decision: keep the unchanged local quadratic as an explicit limited comparator,
not a finite-domain conditional-surface predictor by assertion. The native
sampling is necessary; choosing Hessian eigenaxes and adding formal Schur
elimination does not itself solve finite-displacement accuracy. The next test
must use audited branches and reserved independent points, not refit these
seen observations and call the fit a prediction.

## What the disagreement does not establish

Nonlinearity, finite-difference effects, branch changes and residual relaxation
can contribute. These observations do not uniquely separate those sources.
Neither differences between the two probe steps nor these training residuals
are physical uncertainty bounds. A partial gradient norm at a candidate is
not proof of stable released curvature, smooth continuation or a TS. These
pointwise energy errors are not measurements of whole-path barrier errors.
No H1/H2 independent success or failure is established here.
This does not invalidate the measured bounded optimizer-initialization result:
BFGS updates its matrix, whereas this comparison deliberately keeps a fixed
original quadratic as an energy model. Acceleration and finite-domain model
accuracy are separate tests.

The public quadratic API is used with an explicitly algebraic zero floor only
to construct the raw mathematical comparator. That is not a measured physical
stability floor. No negative eigenvalue is clipped or pseudoinverted and no
post-relaxation released block is certified. Original ABACUS/PBE/100 Ry/full
10-au DZP/Gamma2x2x2 inputs, ordinary NEB0.10, live runtime and jobs are unchanged.

## Reproduce

From the repository root:

```text
python benchmarks/hfo2_channels/20261008/G3_quadratic_agreement_E165_20261015/analyze.py --output a-new-output.json
python -m pytest -q tests/test_hfo2_G3_quadratic_agreement.py
```

Existing output files are refused. Original curvature/canary inputs and both
native prefix reports are hash-pinned. No raw SCF reparse, private author array
or physical +0.005 holdout is accessed. See [complete original curvature](../G3_complete_curvature_E159_20261014/README.md),
[native candidate prefixes](../G3_face_prefix_audit_E163_20261015/README.md) and
[atomic/cell signatures](../G3_atomic_cell_signature_E164_20261015/README.md).

## Next discriminating control, not yet submitted

Within each existing face's remaining training allowance, prepare four same-Q
frozen-R static controls at (+/-0.1,0) and (0,+/-0.1). R=0 means the original
nonstationary reference's omitted coordinates; Q may itself contain open-cell
motion, so this is not a fully fixed-cell experiment. The original Q0/R0 SCF
is already cached. Matched released candidates are already audited. The purpose
is to distinguish disagreement of the frozen finite-Q quadratic from the
effect of the declared relaxation protocol, without claiming global minima.

If admitted, these consume four of each face's 56 uncommitted training
structure/node slots, leaving 52, with all eight independent slots untouched.
They are training diagnostics selected after seeing these training residuals,
not independent tests. No new Hessian center, material, physical strain,
electronic parameter or complete NEB chain is added. Do not alter the immutable
E161 policy files to backdate this later reservation. The same original six
physical input hashes and fixed mechanical boundary apply.

Admission must wait for an existing allocation to finish and an explicit
resource/budget check; at most two committed study allocations remain allowed.
No new frozen-control geometry, input deck, job or reservation is created here.
Signed frozen/released energy differences must be retained even if unfavorable;
do not force a candidate onto a different branch or promote a last trial.

## Verification

307 related checks pass locally and 307 in an independent archive, including
fourteen new case checks and the existing quadratic API regressions. Import
roots are verified, and the independent CLI reproduces the entire saved
analysis byte-for-byte (SHA c2e51c37b185a00676a6ff21660ea132384c0a22cfac2cebd18e30a5e5a9b745).
The frozen source/native-report byte checks are not replaced by the arithmetic
replay tolerances. Tested tree fa0fa344da8fa598a1493fdb4f0b71795466097d and exact
source/archive/test evidence are in validation.json. This is not a new
whole-repository regression, raw SCF re-audit or independent material forecast.
