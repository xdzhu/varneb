# Pairing a stationary well and bottleneck response

`vcneb.stationary_gap.restricted_stationary_gap_response` pairs the
[restricted stationary branches](CONTROLLED_STATIONARY_BRANCH.md) at the
**same prescribed physical parameter**. It is standard local harmonic response,
not a new theorem or a completed B2--B5 material predictor. No calculator is
called, no branch is selected and no activation barrier is certified.

## Shared physical declarations, distinct local charts

Each `ControlledStationaryModel` contains a validated quadratic branch,
one `StationaryResponseContract`, a positive external scale, explicit parameter
bounds, a displacement radius and its source audit hash. The common contract
binds input hashes, mechanical family, fixed boundary definition, control name
and unit, anchor condition, energy zero, formula normalization and pressure.
The fixed-boundary hash must include the substrate reference, open freedom,
loads/field and other invariant settings; do not hash a changing strained cell
and then call it an independent boundary change.

Both contracts must match. The initial branch must have represented internal
index 0 and the bottleneck index 1. Raw E supports P=0; at finite pressure both
models must use same-reference H=E+PV, including its derivatives. Hashes and
these declarations do **not** independently audit the supplied physical data.
HfO2 production still uses its original P=0/E=0/100Ry/full10-auDZP contract.

Local centres can have different cells, atom/strain charts and positive
parameter scales. For centre a, `t_a=L_a*(x-x0)`; evaluate both models at the
same x, not the same numerical t. Their Hessians and physical gradients must
each be in their own correctly declared Angstrom chart. Reference mapping,
strain normalization and probe construction remain externally audited inputs.

## Preserve offsets; do not subtract large totals to obtain response

Let e_a(t) be the branch energy change relative to its raw centre energy,
including the nonstationary reference correction. The caller supplies the
independently recorded raw centre gap G=E_S,ref-E_I,ref. Its agreement with
the paired raw totals is checked with an explicit identity tolerance, not
treated as prospective accuracy or a changed SCF threshold.

```text
D(x0) = [G + e_S(0) - e_I(0)] / n_formula_units
d_a(x) = e_a'(0)*t_a + curvature_a*t_a**2/2
delta_D(x) = [d_S(x) - d_I(x)] / n_formula_units
D(x) = D(x0) + delta_D(x)
dD/dx = [L_S*e_S'(t_S) - L_I*e_I'(t_I)] / n_formula_units
d2D/dx2 = [L_S**2*curvature_S - L_I**2*curvature_I] / n_formula_units
```

The output separates G, the stationary-anchor correction, D(x0), each energy
increment, delta_D, the signed D(x) and both derivatives. Increments are
computed from small local coefficients rather than subtracting large total
energies. A common raw-energy shift therefore does not change the response.
Negative gaps remain negative; they are not clipped into apparent barriers.

For the same **actual** stationary IS and connecting TS, D is the activation
energy and these derivatives describe its branch response. A restricted
quadratic solution alone cannot establish those premises. In particular,
adding delta_D to an audited ordinary-NEB barrier is a separate calibration:
record that choice and check the stationary-anchor mismatch instead of hiding
the offset or silently replacing the audited barrier. Unconverged G1 peaks
are not certified TS centres for this purpose.

## Domains and failure are explicit

Both the stationary anchor and requested x must be in each declared parameter
interval and displacement radius. Out-of-domain, inconsistent-reference,
wrong-index and nonfinite responses are rejected; no extrapolation repair,
recentering, pseudoinverse or physical-input retuning occurs. These bounds are
caller-supplied necessary checks, **not** proof that the response remains in
the actual probe convex hull or that a quadratic is accurate.

The result keeps both complete `StationaryQuadraticPoint` diagnostics:
represented residual, unrepresented but physically movable gradient, and
mechanically clamped reaction. Frozen B2-type omissions can have nonzero
movable gradients; do not erase them or call projected stationarity full
stationarity. A clamped reaction is not a failed open traction. A negative
curvature outside measured/retained directions is not ruled out by this API.
Independent DFT, anharmonicity, connection/continuity, alternate bottlenecks,
branch changes, sampling error and lowest-channel coverage still need tests.

## Reproduction and current limits

Run from the repository root with a new output file:

```sh
python -m scripts.check_stationary_gap --output /new/paired-check.json
python -m pytest -q tests/test_stationary_gap.py tests/test_stationary_branch.py tests/test_biaxial_curvature.py
```

Eight analytic combinations of frozen/stable release, rotated charts and
equal/distinct external scales test 24 points against separately assembled
stationary solves and finite gap derivatives. Tests also cover common energy
shifts, actual biaxial geometry with distinct centre cells/scales, source
identity, bounds, ownership, omission diagnostics and signed negative gaps.
These are implementation residuals, not DFT uncertainty or material accuracy.
The receipt and clean/HF replay are kept in
[the E032 delivery](../benchmarks/hfo2_channels/20261008/stationary_gap/README.md).

The registered [prediction controls](PREDICTION_CONTROLS.md), same-boundary
G2 training, training-only basis/branch choice, complete B2--B5 freezing and
independent holdouts remain required. This module does not authorize those
jobs or claim a favorable HfO2 selectivity window, quantitative model advantage
or JCTC-ready manuscript. It is the paired-response step, not the end state.
