# Sampled physical-gradient work along a declared path

`vcneb.path_work.sampled_path_work` is calculator-independent bookkeeping of
physical covectors, not another optimizer or energy decomposition theory.
Use a conjugate coordinate/covector pair and raw E (or consistently H=E+PV).
The routine neither infers that chart nor certifies a calculator or uncertainty.

For each segment, delta x = x(i+1)-x(i) and gbar = [g(i+1)+g(i)]/2. With
complete mutually orthogonal subspaces U(k), the signed discrete work is

```text
W(i)    = gbar . delta x
W(k,i)  = (U(k).T gbar) . (U(k).T delta x)
sum_k W(k,i) = W(i)
R(i) = [E(i+1)-E(i)] - W(i)
```

Each subspace is invariant to rotations of its basis, but choosing a different
span or chart changes its signed work. Work can be negative. It is not a
state-function energy associated with one mode, a coupling coefficient or a
prediction of a path after orthogonal relaxation. Do not integrate NEB
perpendicular/spring forces in place of the physical energy gradient.

Keep mode-mass coordinates and scaled-cell coordinates as separate blocks:
their units may differ, but their covectors must be dual and no basis may mix
incommensurate units. Atomic gauge translations are another explicit block.
Complete/rank-preserving partition checks cannot establish local stability,
conditional minimization, a saddle index or numerical DFT accuracy.

R(i) and its cumulative sum are **diagnostics**, not physical error bounds.
Even exact physical gradients at coarse endpoints can miss nonlinear variation
between them. Small final closure can result from cancelling interior residuals.
No force threshold or quadrature residual is promoted to barrier uncertainty.

The [actual four-band HfO2 example](../benchmarks/hfo2_channels/20261008/clamped_path_work_E092_20261011/README.md)
reparses36portable terminal SCFs in the original100Ry/full10auDZP contract.
Its atom/cell pullback is independently checked against native Cartesian
force/ASE-stress contractions. All original T/Cmma spans, ranks and reference
registrations remain fixed; unknown/held-out receipts reject before access.
No calculator, scheduler write or extra physical Hessian is involved.

```text
python -m pytest -q tests/test_path_work.py tests/test_hfo2_clamped_path_work.py
```
