# GaN B4→B1: a real two-dimensional mode-resolved enthalpy study

Status: active experiment, **not yet a computed two-dimensional surface**.
Pressure: 45.7 GPa. The four-atom Ga2N2 cell contains two GaN formula units.
Every plotted energy must be the enthalpy `H=E+PV` at this pressure and must
state whether the reference is B4, B1, or the near-TS center. **All new GaN
mode, surface and path calculations use the original VASP production contract:
PBE, Ga_d+N PAW, ENCUT=600 eV, Gamma 8×8×6, ISYM=-1, SYMPREC=1e-4,
EDIFF=1e-7.** Keep the original POTCAR and k-point convention, atom order,
pressure and VASP version. Do not change ENCUT for a diagnostic, a pilot or an
individual failed point, and do not mix 600-eV and 1000-eV enthalpies or mode
vectors as though they were one quantitative surface.

## Why this is still required

The existing five-backend figure is five one-dimensional, 29-total-image
paths. A separate **1000-eV historical diagnostic** found a near-stationary,
locally index-one joint atom–strain candidate (`|g|=0.00156 eV/Å`, lowest two
curvatures about `-4.03,+2.56 eV/Å²`) and two-basin links. It does **not**
compute a two-dimensional surface, and it cannot supply the quantitative TS
mode basis or saddle certification for the 600-eV production result.

An offline, hash-checked **same-600-eV geometry diagnostic**, not a computed
enthalpy surface or TS certificate, compares the final 600-eV chain with its
image-15 600-eV joint Hessian in
`benchmarks/numerical_integrity/gan_2d_joint_mode_feasibility_20260927.json`.
The fixed plane spanned by the image-15 negative mode and one path-fitted
transverse joint direction misses the chain by up to **0.290 Å** (95.85% of
the squared displacement norm captured). The image-15 lowest two eigenmodes
miss by **0.499 Å**; even the unconstrained best rank-2 joint plane misses
by **0.096 Å**. The periodic gauge margin is 0.387 and cell-rotation residual
only `1.5e-6`. The image-15 full-gradient norm is `0.07093 eV/Å`, so its
negative Hessian mode must **not** be called the certified TS unstable mode.
The geometric gate rejects a single fixed two-mode plane for the whole path;
it supports testing path-adapted coordinates, not claiming the tube is
computed already.
The existing 600-eV finite-difference audit also reports a maximum
energy-gradient/stress-gradient discrepancy of `0.02372 eV/Å` at a 0.02 Å
step (and `0.02320 eV/Å` at 0.01 Å). This is a numerical consistency limit
that must be reported and diagnosed **without changing ENCUT** before making
sub-0.01 eV/Å saddle-gradient claims.

## Two complementary computed surfaces, with distinct claims

1. **Local TS joint-mode cut.** Starting from the 600-eV production path's
   barrier neighborhood, refine the saddle and compute its joint atom–strain
   Hessian **entirely at 600 eV**; remove translations and choose the negative
   eigenvector `u` and a declared orthogonal transverse eigenvector `v`.
   First verify the stationary gradient, Hessian index and two-basin links
   under this same contract. Define a *frozen* local slice
   `H(q_u,q_v)=H[Q_TS+q_u u+q_v v]`. Only after these gates pass, preflight a
   small 3×3 static pilot at `q_u,q_v∈{-0.10,0,+0.10} Å`, reusing the center
   only if its input and output are genuinely identical. The plot is a
   local *enthalpy cut*, not a full-path PES or a constrained minimum. Check
   center/axis energy–gradient consistency, cell positivity, minimum
   distances, SCF termination, input hashes, and the expected negative/
   positive curvatures before any interpolation. Only then add held-out
   interior points where the pilot fails a predeclared meV/GaN error target.
2. **Whole-path path-adapted tube.** Let `s` be the audited VCNEB joint-metric
   arc length and let `v_perp(s)` be an explicitly transported, normalized
   joint atom–strain transverse direction. The actual 600-eV VCNEB images
   define `q_perp=0` by construction; off-path statics at controlled signed
   `q_perp` produce `H(s,q_perp)` under the **same 600-eV production contract**.
   Start with nine preselected anchors (including both endpoints and the
   barrier neighborhood), two off-path points per anchor, and same-protocol
   center values; adapt only cells with large held-out interpolation error
   or branch changes. A path-tube chart is *not* a global two-phonon plane:
   `s` is a collective reaction coordinate, and the transported transverse
   direction can change with `s`. Plot the mode/strain amplitudes and
   reconstruction residual separately along all 29 images.

The two cuts answer different questions and should be separate panels or
figures. Their **same electronic contract** does not make their coordinate
charts interchangeable: the local TS axes and the path-adapted axes must
remain labeled separately. An optional conditional relaxation at fixed `(s,q_perp)`
is a later, explicitly labeled calculation: frozen, relaxed-orthogonal, and
strict-subspace barriers remain different quantities.

## Gates and stop conditions

- Before any new DFT submission, verify every generated cell and the
  calculation-setting hashes (INCAR/KPOINTS/POTCAR) against the original
  600-eV path; POSCAR is expected to differ across sampling points. The
  aborted 1000-eV array
  `27791808` is excluded: two cases completed, six failed with a VASP
  direct/reciprocal Bravais-classification conflict; no result is used.
  Diagnose that input-geometry issue and demonstrate a same-contract
  preflight on the actual proposed cells before launching a new surface.
  Prior equivalent-basis and explicit-identical-k-mesh canaries did not
  resolve related conflicts; do not re-label those as validated repairs.
- First pilot: at most eight independent TS statics after the 600-eV saddle
  and mode basis are audited. If any source/input hash,
  geometry, SCF, force/stress, or enthalpy-gradient check fails, diagnose the
  cause before adding grid points. Do not change electronic settings per
  point to rescue a plot.
- Local figure: show every measured DFT point; shade/interpolate only inside
  their convex hull, report independent holdout error and do not infer a
  barrier from color contours alone. The true candidate is checked by a
  same-protocol gradient and one-negative-direction test.
- Full path: `q_perp=0` must reproduce the archived 600-eV path energies;
  endpoint, atom order, gauge and strain coordinates remain explicit. If
  transverse branch identity changes, split the chart instead of smoothing
  over a discontinuity. No global fixed-two-mode claim unless its off-plane
  residual is independently reduced below the predeclared 0.05 Å gate.
- Completion for the requested mode analysis means: a raw-audited TS local
  2D cut, a raw-audited whole-path 2D tube, all-image mode/strain/reconstruction
  plots, a same-protocol comparison to the original 1D enthalpy barrier, and
  figure source data. HfO2 and bilayer hBN come **after** this GaN/BTO closure.

For a stationary enthalpy minimum or an index-one saddle, the derivative
`dH/ds` along a smooth MEP is zero at that point. An index-one saddle has
negative reaction-direction curvature; stable endpoints have positive
curvature. The GaN barrier is therefore not two independent barrierless
transitions joined at a stable intermediate: the shared TS is unstable.
This statement concerns the 0 K enthalpy path, not a finite-temperature
free-energy profile or the slope of a discrete image-number plot.
