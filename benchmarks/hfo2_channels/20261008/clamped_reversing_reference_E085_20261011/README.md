# E085: completed reversing-band fixed-reference and external-work analysis

Zero new DFT. This applies the existing T chart, ranks3/4 and all four E023
Cmma registrations to the ordinary-converged reversing band28692776,
step5/.09218242631 eV/Angstrom. Original reference source files and geometric
frame/order/lift are unchanged. This is training description, not a new
phonon calculation, G3 mode/bottleneck selection or prospective forecast.

At the sampled peak4, using the same translation-free atomic chart:

| Reference span | Rank | Squared displacement (%) | Squared gradient (%) |
|---|---:|---:|---:|
| T geometric triplet |3|26.2915|0.0320|
| Two lowest complete T optical doublets |4|16.2196|10.0215|
| All four registered Cmma frames |4|46.5637--46.5639|65.1046 at displayed precision|

The common internal gradient norm is0.06109916898 eV/(sqrt(amu)*Angstrom),
Cmma omitted norm0.03609271285--0.03609271465. Displacement and gradient
fractions have different denominators; neither partitions energy or measures
released-model/force/barrier accuracy. Atomic and open-cell gradients remain
separate. Source LDA directions are not new PBE local modes/curvatures.

The frozen, fixed-open-row substrate work at peak4 is-27.183577 meV/HfO2
per unit dimensionless strain; relative to PO it is201.230946. It is a
declared-coordinate partial derivative, **not a stationary barrier slope**.
The almost untilted geometry gives tiny affine/fixed-open-row work difference
(1.4653e-11 eV/cell/unit strain), not evidence of full TS stationarity.

## Evidence and reproduction

- `displacement.json`: unchanged reference analyzer; nine portable native
  SCF frames reparsed, prescribed inverse proper cyclic frame on copies.
- `gradient.json`: dual physical covector transport, all reference frames and
  norm/translation diagnostics; no force/source-mode ASR correction.
- `external_work.json`: raw tensile-positive ASE stress and exact original
  clamped deformation Jacobian; no finite-strain/Hessian calculation.
- `shared_initial.json`: three completed training bands' cached PO native
  INPUT/KPT/STRU/logs and EFS compared, ordered geometric serialization
  differences measured explicitly; no repair or recomputation.
- Actual receipt632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41
  pins E083's historical HF42-SCF/six-file/338-runtime audit. Local analysis
  reparses portable data, not undistributed orbitals or every earlier SCF.

Use the two original commands in E071/E080 with the E083 case/receipt and
fresh output paths, keeping original seed/reference/source arguments. The
six pinned author files remain outside Git; raw or transformed eigenvector
arrays are not redistributed. Frozen-work reproduction is:

```text
python -m scripts.analyze_hfo2_clamped_external_work --case benchmarks/hfo2_channels/20261008/clamped_reversing_terminal_E083_20261011 --audit-sha256 632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41 --output NEW_REPORT.json
```

The shared-work helper now admits this actual pinned terminal receipt; it
still rejects arbitrary receipts before reading labels and does not admit
the+.005 test condition. `validation.json` records executed tests, independent
archive replay, the actual manuscript build and inspected layout.

An initial exact floating-array comparison of the shared geometry failed on
the previously known4.44e-16 Angstrom serialization difference. Exact native
SCF hash/EFS comparisons remain; geometry is checked and its difference
retained rather than modified. The failed JUnit receipt is preserved alongside
the corrected tests. This is not a material or SCF failure.

The manuscript adds three descriptive Table6 rows and their absolute-gradient
context. No sampled peak is promoted to a continuous saddle, no reference
is dismissed from its projection alone, and no independent model advantage
is asserted. Full G2/G3/strong-control/frozen-holdout evidence remains missing.
