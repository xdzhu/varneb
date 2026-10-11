# E081: clamped mode-calculator flow, not a material conditional surface

The old `CalculatorModeEvaluator` accepted only a six-symmetric-strain
reference chart. That chart cannot represent the exact G2 boundary: two
prescribed substrate vectors, with the third cell row released (including
tilts). Replacing the boundary would change the scientific experiment.

`ClampedModeCoordinates` therefore reuses the existing
`ActiveJointCurvatureCoordinates.for_clamped_plane` geometry and physical
gradient, without introducing another backend or changing production inputs.
The same ASE calculator factory and native SCF validator can now feed fixed-Q
conditional optimization. The old chart and serialized contracts remain
unchanged. See `docs/CLAMPED_MODE_EVALUATIONS.md` for units and use.

## Evidence and limitations

`implementation_verification.json` records an actual local verifier run:

- Eight oblique/rotated, tilt/no-tilt, zero/nonzero-pressure EMT cases use
  independent central energy differences: 232 EMT evaluations in total.
  Maximum work-gradient disagreement is below 1e-9 eV/Angstrom.
- Exact-coordinate cache reuse and serialized recovery are exercised.
- E079 PO-to-T and E068 PO-to-M provide 18 distinct portable native terminal
  frames. Each ordered cell/position reconstructs in one actual T-endpoint
  reference within 1e-12 Angstrom. Energy and active-gradient replay exactly
  match the unchanged primitive. This is not an independent native DFT
  energy/force finite-difference measurement.
- Native six-asset/runtime audits were completed in E068/E079, not repeated
  by this local verifier. No new DFT, production deployment, scheduler write,
  material Hessian, mode selection, surface or holdout access occurs.

The tests also exercise conditional safeguarded BFGS with explicit fixed Q,
translation gauge and pre-calculator geometry rejection; wrong substrate,
backend geometry mutation, malformed/asymmetric stress and changed resume
contracts do not become valid results. The analytic EMT gradient tolerance
does not replace production NEB 0.10 eV/Angstrom or endpoint settings.

The first injected-asymmetric-stress test exposed ASE's raw-to-Voigt tensor
averaging. The new adapter checks raw stress before that conversion; it does
not symmetrize backend output. The first conditional integration test omitted
the already mandatory geometry validator; the test was corrected without
changing the optimizer or material settings.

`local_tests.xml` records 70 final related tests (including 25 new adapter
tests); additional full-project and actual clean-archive verification, source
identity and direct process exits are recorded only when complete in
`validation.json`. Static `units_audit.json` has zero heuristic findings;
it is not an uncertainty bound or proof of DFT convergence.

## Reproduce without launching DFT

From the repository root, use fresh output paths:

```text
python -m scripts.verify_clamped_mode_evaluator --workdir NEW_WORKDIR --output NEW_REPORT.json
python -m pytest -q tests/test_clamped_mode_evaluator.py tests/test_mode_evaluator.py tests/test_clamped_joint_curvature.py tests/test_mode_surface.py
```

G2 training remains incomplete. G3 selection is still deferred until the
complete registered training set, and the reserved +0.005 control remains
unread. The original 100 Ry/full 10 au DZP/PBE/P0/Gamma 2x2x2 six-file
contract, ordinary NEB threshold and two-allocation limit are unchanged.
The whole-grid interpolation gate is also unchanged: eight independent G3
points must not be described as validating all 64 cells of a 9x9 grid.
This is a tested implementation milestone, not H1/H2, predictive advantage
or a completed JCTC manuscript.
