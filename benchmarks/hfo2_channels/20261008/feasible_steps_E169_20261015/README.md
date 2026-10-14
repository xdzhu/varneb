# E169 — an opt-in ball-feasible step experiment, not material validation

Previous E168 was PROGRESS, pushed as0d65557b: the real preserving point13
replay identified35 radius refusals and a0.0723-eV/Angstrom remaining tangent
gradient. This milestone implements a specific next hypothesis without
changing the old solver, Hamiltonian, target, live runtime or domain.

The new independent module uses the existing ModePlane, metric-normalized
released directions, fixed-control checks and ConditionalModePoint. It takes
one explicitly audited start. Proposed steps are projected into the declared
released ball; non-descent projected quasi-Newton steps fall back to projected
steepest descent, with guarded Armijo and standard positive-curvature BFGS.
The full released gradient, not boundary KKT stationarity, controls success.
No backend-specific code, physical input, new default or source repair is used.

See the [API contract](../../../../docs/FEASIBLE_CONDITIONAL_RELAXATION.md).
In particular a released radius is not silently identified with a total
displacement limit. The original guard must remain in the material experiment.

## Frozen analytic cohort

Before execution, `benchmark.py` defines four cases, two arms, one common
radius0.5, norm target0.03, trust radius0.05 and unique-call cap12. It changes
no coefficient after seeing the comparison. Both arms use the same exact-
coordinate analytic cache and seed. All four cases remain in the result.

| Analytic case | Original calls | Feasible calls | Original guard refusals | Feasible guard refusals | Outcome |
|---|---:|---:|---:|---:|---|
| Interior quadratic | 4 | 4 | 0 | 0 | Both reach target |
| Near-wall, misaligned positive seed | 7 | 6 | 7 | 0 | Both reach target |
| Nonlinear valley | 9 | 9 | 0 | 0 | Both reach target |
| Exterior quadratic minimum | 11 | 8 | 42 | 0 | Neither reaches full-gradient target |

These counts are unique analytic oracle evaluations, not DFT. Callback counts
and complete trial gradient norms are retained in `analysis.json`. The
exterior case has no candidate in either arm; its failed costs are not a
speedup ratio. The two no-improvement controls are not dropped, and the one-
call improvement on the artificial near-wall case does not establish material
benefit or a new optimization theory. A separate analytic ASE E/F/stress
adapter test verifies that the conservative factory budget still propagates.

## Reproduction and next gate

```text
python benchmarks/hfo2_channels/20261008/feasible_steps_E169_20261015/benchmark.py --output a-new-output.json
python -m pytest -q tests/test_feasible_relaxation.py tests/test_feasible_steps_benchmark.py
```

Actual HfO2 testing still requires a prospectively frozen matched start,
unchanged full/released radius identity, force target, Hessian seed, original
six physical hashes and explicit unique-call/reservation accounting. A new
runtime must be separate from the live E161 source; its old provenance pins
must not be bypassed. The already prepared E166 frozen statics keep priority.
No new material continuation, geometry, independent label, Hessian center,
Slurm task or training-slot reservation is created by this milestone.
The full stable-face/strong-control independent-prediction goal remains open.

## Verification

The final ten-module suite passes218 checks locally and218 from a fresh
independent Git archive, with zero failures/errors/skips;19 checks are new.
These cover fixed Q, nonunit metric/frozen controls, seed validity, domain
refusal before evaluation, full-gradient rejection of boundary stationarity,
budget exception propagation, an analytic ASE E/F/stress factory budget and
the complete four-case control cohort. This is not a whole-repository or
native material-validation result.

The independent CLI reproduces the complete saved benchmark byte-exact:
SHA20c5d743f267be926bb2382f30aeec32b2a46172cb74f74e9223c9ff76153a1d.
Tested tree4b1a61ca7459eee56a3415e50a0e724502afa005 and source/archive/test
hashes are retained in `validation.json`. Existing mode_surface/default code
and all19 normalized scientific code pins from E166 remain unchanged.
No live source is overwritten; no material job, extra training reservation,
physical +0.005 structure or label is created. The prospective DFT comparison
and the full independent-prediction goal are still unfinished.
