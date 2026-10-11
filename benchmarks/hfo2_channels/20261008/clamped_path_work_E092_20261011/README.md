# E092: signed native path-work audit, not mode energy allocation

Use only the four currently admitted ordinary-converged G2 training bands:
zero-strain PO-to-M, PO-to-T and reversing flip, and +1% PO-to-T. Do not read
the ongoing preserving/+1% M chains or reserved +0.5% holdout for this analysis.
Retain the same T mass/translation gauge, rank3 geometric triplet, rank4 two
complete optical doublets and all four rank4 Cmma frames103/106/114/127.
Author mode arrays remain outside the repository; their original pins are
verified by the existing reference audit. No author code or calculator runs.

The numeric helper integrates physical energy covectors against successive
coordinate differences with endpoint trapezoids. Each signed contribution
is an explicit path-dependent work in a declared subspace; it is **not** a
state-function energy belonging to that mode. Atomic retained, omitted and
translation work plus the original permitted open-cell work reproduce the
full discrete contraction. Basis rotations within a span change no result.

Raw native Cartesian forces/stress contracted directly with the native
fractional/cell polygon provide an independent cross-check of the dual
common-T atomic chart and native open-cell chart. No spring or NEB forces are
integrated. Each band uses its own matched initial PO+ cell and raw P0 energy.
The complete four-band work profile is not a complete eight-band G2 network.

Energy-minus-work discrepancies are retained per segment and cumulatively.
Coarse sampling/nonlinearity and DFT force/stress error are not separated by
this diagnostic. A discrepancy is not a barrier error bound, failure of a
mode model, proof of a continuous TS or justification for changing inputs.
The ordinary0.10 criterion and 100Ry/full10auDZP contract stay unchanged.

Reproduce using the same source/root/registration/variants/Gamma/seed arguments
as `scripts.analyze_hfo2_clamped_reference`, but with the module
`scripts.analyze_hfo2_clamped_path_work` and a fresh output path. Exact terminal
audit pins are in the four reports. Unknown receipts reject before source
access. Public author-file provenance is in `../cmma_reference/README.md`.

```text
python -m pytest -q tests/test_path_work.py tests/test_hfo2_clamped_path_work.py tests/test_hfo2_clamped_mode_gradient.py tests/test_clamped_mode_evaluator.py
```

Synthetic quadratic/rotation/refinement checks are implementation checks,
not HfO2 predictive accuracy. The four-band test reparses36portable native
terminal SCFs, not all earlier undistributed SCFs/PP/orbital/runtime files.
These36frame records represent34distinct raw log/STRU identities because
three zero-strain paths reuse one initial PO+ calculation. They are not36new
DFT calls or independent statistical replicates.
G3/reference/model selection, extra Hessians, complete forecast freeze and
independent prediction evidence remain pending. This analysis makes no
physical/scheduler mutation and does not change production workers.

## Actual measured diagnostic

All values below are meV/HfO2 from PO+ to the highest sampled image. These
are signed trapezoidal work and raw sampled energy differences, not separately
assigned mode energies or fitted barrier values.

| Completed training band | Raw sampled peak | Full work | Energy minus work | T-triplet work | Omitted atomic work | Open-cell work |
|---|---:|---:|---:|---:|---:|---:|
| zero strain PO-to-M |98.458122|77.560379|20.897743|-0.111271|68.756917|8.914733|
| zero strain PO-to-T |84.416624|75.955937|8.460687|2.666921|72.789643|0.499373|
| zero strain reversing flip |338.968965|303.693643|35.275323|195.447066|105.539695|2.706882|
| +1% strain PO-to-T |96.763650|87.345393|9.418257|1.367408|85.773095|0.204891|

Rigid-translation work is also retained in each report (at most3.90e-9 in
these units at the displayed peaks), not silently repaired as source ASR.
Independent native versus dual-chart work agrees within1.12e-16eV/cell.
Endpoint energy-minus-work is only3.006016/0.728018/0.115373/1.997185 for
these four bands: endpoint closure can hide much larger interior discrepancies.

For the **same** sampled path and gradient information, Cmma retained atomic
work ranges over all four registrations are55.618955--55.619265,
62.066996--62.068212,63.846564--63.846572 and67.984207--67.985028 respectively.
The two complete T optical doublets give12.359091/11.900861/0.782574/11.007299.
The two rank4 controls and rank3 T triplet remain distinct, with all omitted
work included. No reference is selected based on these results and no
conditional-branch or prediction-accuracy ranking is inferred from them.

Thus high displacement coverage or final endpoint closure cannot substitute
for an energy-resolved validation through the bottleneck region. This supports
retaining the planned measured conditional response and fair strong-reference
comparison; it does not prove a new method wins or mandate changing electronic
parameters. Quadrature error versus physical force/stress error remains
unseparated until appropriate measured sampling evidence exists.

## Final validation and live state

Local focused tests pass67; actual clean full regression passes2187 with2
explicit existing skips,0failure/error and314existing warnings in361.051
JUnit seconds. No HF pytest or software is installed. Each of the four
actual material analyses also executes from the clean archive; all exported
reports are byte-identical. Final independent delivery tests pass67/0failure/
0skip in4.114seconds, actual exit0. Eleven then-present case files and all
four new code/test files match that archive; code/test bytes also match the
full regression archive. Sources and material records do not change afterward.
Skips are the unavailable optional ASE ABACUS reader and undistributed licensed
VASP POTCAR, not passed backend validations. Receipts identify actual runs.

At20:02:30CST the original two study allocations remainRUNNING/32CPU,
stderr0/no failure: +1%M28727368 step2/.525706 and preserving28727382
step2/.140697eV/A. The two +1%flip starters stay held/nullDependencies.
G2 remains4/8. This milestone performs0schedulerwrites/admissionreplays/new
chains or physicalDFT; production workers continue normalDFT. No G3 choice,
extraHessian or+.005physicalholdout access. Main manuscript/PDF and original
five user-dirty tracked files remain unchanged. Full scientific goal staysactive.
