# E027: a controlled substrate parameter, not a released degree of freedom

The [coordinate contract](../../../../docs/BIAXIAL_CONTROL_CURVATURE.md) adds
the external biaxial control to a local sampling chart. At each prescribed
strain, the endpoint/NEB internal mechanical boundary remains unchanged.
This is a derivative implementation prerequisite, not a new elasticity theorem,
verified hafnia prediction, conditional DFT surface or TS certificate.

## Measured checks

- Local isolated full regression: **1073 passed, 2 skipped, 313 warnings**;
  107.29 seconds, baseline commit `407076a` plus the four recorded code/test files.
  The user's dirty/untracked files are not in that snapshot.
- Focused regressions: 130 passed / 6.42 seconds; 28 new controlled-coordinate cases.
- The same byte-checked NumPy/ASE checker on local/HF: eight Cu/EMT gradient
  configurations, 256 evaluations per host, ten real **uncomputed** HfO2 seed
  geometries and two mixed-derivative step sizes. Maximum full-gradient error
  was 1.50e-9 locally and 1.84e-9 eV/A on HF; these are implementation residuals,
  not DFT energy/force uncertainties. Mixed reciprocity defects stayed below
  1.3e-10. Complete metrics and source hashes are in
  [validation_delivery.json](validation_delivery.json).

The HF `icu` environment has no pytest. The attempted pytest run stopped with
`No module named pytest`; nothing was installed. HF therefore ran the explicit
NumPy/ASE check, **not the full pytest suite**. The helper first failed under
direct nested-file execution because the repository root was absent from
`sys.path`. Running its unchanged bytes through stdin from the source root
fixes that invocation, not the coordinate algorithm.

Reproduce the dependency-free checker **from the repository source root**:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python - < benchmarks/hfo2_channels/20261008/biaxial_control/hf_numpy_compatibility_check.py
```

PowerShell equivalent, also from the source root:

```powershell
Get-Content -Raw benchmarks/hfo2_channels/20261008/biaxial_control/hf_numpy_compatibility_check.py | python -
```

The actual HF check used a fresh immutable source bundle under
`/public/home/iai806/abacus/agent-runs/20261009-varneb-biaxial-control-e027/source`
and base64-decoded stdin to verify the checker bytes explicitly. It did not
call ABACUS, change production sources, install packages or submit jobs.
The Cu pressure in the checker is for derivative validation only; production
HfO2 remains P=0 with the registered 100 Ry / full 10-au DZP contract.

## What remains

The augmented 40-coordinate HfO2 chart contains 39 internal coordinates plus one
prescribed input. The latter must not relax or be counted as a physical unstable
mode in a fixed-strain TS test. Internal release remains stability-gated; sampled
couplings do not establish unsampled positive curvature. No B2--B5 forecast batch,
material Hessian, independent barrier prediction or favorable H1 window is
completed by these checks. G1 -> G2 -> G3, the finite budgets and the held-out
condition are unchanged.
