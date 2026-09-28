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
At the 0.02 Å step the three diagonal-strain discrepancies are nearly
isotropic (about `0.022–0.024 eV/Å`), whereas the largest atomic-coordinate
discrepancy is `0.00674 eV/Å`. This pattern is compatible with a stress-side
basis-set/Pulay contribution, but does **not** prove that diagnosis; VASP's
[Pulay-stress documentation](https://vasp.at/wiki/Pulay_stress) explains why
cell derivatives are especially sensitive to a finite plane-wave basis.
For this study the response is to carry a same-600-eV numerical error budget,
not to retune the cutoff or combine unequal-cutoff energies.

The current CPC draft's GaN joint-mode panel still depicts the segregated
1000-eV diagnostic. It is **not** the final same-protocol figure and must be
replaced or removed before submission; the historical calculation is retained
only for provenance. Endpoint Γ modes, local joint modes and the 1D path must
all be compared under the 600-eV contract in the revised analysis.

## Two complementary computed surfaces, with distinct claims

1. **Local TS joint-mode cut.** Starting from the 600-eV production path's
   barrier neighborhood, refine the saddle and compute its joint atom–strain
   Hessian **entirely at 600 eV**; remove translations and choose the negative
   eigenvector `u` and a declared orthogonal transverse eigenvector `v`.
   Verify the stationary gradient, Hessian index and two-basin links
   under this same contract before claiming a certified TS. Define a *frozen* local slice
   `H(q_u,q_v)=H[Q_TS+q_u u+q_v v]`. Only after these gates pass, preflight a
   small 3×3 static pilot within geometry-safe displacements, reusing the center
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

As of 2026-09-28, the 600-eV local candidate has one negative joint Hessian
direction at a 0.02 Å full finite-difference step and a 3×3 frozen pilot at
`q_u=±0.02 Å`, `q_v=±0.0125 Å`, plus four independent axial half-steps. The
half-step quadratic-model error is below `0.05 meV/GaN`, but the full
energy–gradient mismatch is about `0.0232 eV/Å` and two-basin links have not
been demonstrated. Thus this is a measured local cut, not a certified TS.
The old 1000-eV coupling figure has been replaced in the CPC draft by a
same-600-eV Hessian/endpoint-Γ panel; the 1000-eV diagnostics remain historical
only and cannot be mixed into production-path claims.

The already completed 600-eV **local** pilot now has a separate plotted
quadratic cut, `paper/VARNEB_CPC/figures/gan_600eV_local_joint_cut.*`.
Within the measured `q_u=±0.020 Å`, `q_v=±0.0125 Å` rectangle, it overlays
eight raw grid statics, four axial half-step holdouts, and the center on the
Hessian-based contour. The maximum grid/model residual is
`0.09054 meV/GaN`; the independent axial half-step maximum is
`0.04974 meV/GaN`, below its predeclared `0.20` gate. A second panel compares
all 12 off-center DFT points with model predictions. This closes the
**local frozen two-coordinate figure and axial validation**, not an off-axis
interior holdout, a relaxed conditional plane, a strict TS certificate, or
the requested whole-path 2D landscape. The source CSV and QA JSON record
the all-600-eV input provenance and the `0.0232 eV/Å` energy–gradient
mismatch. Editable PDF/SVG and PNG are tracked; the 600-dpi TIFF is
reproducibly generated locally but repository-ignored.

The two cuts answer different questions and should be separate panels or
figures. Their **same electronic contract** does not make their coordinate
charts interchangeable: the local TS axes and the path-adapted axes must
remain labeled separately. An optional conditional relaxation at fixed `(s,q_perp)`
is a later, explicitly labeled calculation: frozen, relaxed-orthogonal, and
strict-subspace barriers remain different quantities.

### Same-600-eV two-sided basin-link pilot (raw-audited 2026-09-28)

The audited local cut contains two static seeds at `q_u=−0.02/+0.02 Å`,
`q_v=0`, lower than the near-stationary center by `0.324827/0.490786 meV/GaN`.
Their raw OUTCAR and original VASP 600-eV input hashes are checked
before any continuation. A bounded **native VASP** `ISIF=3`, `IBRION=2`,
`NSW=10` pilot has been staged for these two seeds only, with the *same*
600-eV/PBE/Ga_d+N/Γ8×8×6/EDIFF/ISYM/SYMPREC contract and
`PSTRESS=457 kbar` (45.7 GPa). The job is Slurm array `27793196` on
`hfacnormal01`, 32 MPI tasks per branch, at most two simultaneous. The
remote work root is
`/public/home/iai806/abacus/agent-runs/gan_600eV_ts_basin_pilot_20260928`;
its submitted manifest SHA-256 is
`7fc24267fbed98edcef5c28de757907745adaf651e99f8481e5b011276d53d63`.
A portable, newline-normalized content copy is
`benchmarks/numerical_integrity/gan_600eV_ts_basin_pilot_inputs_20260928.json`.
Both array elements completed with exit code zero (12:33 and 10:37 elapsed).
The raw-output audit is
`benchmarks/numerical_integrity/gan_600eV_ts_basin_pilot_audit_20260928.json`:
all ten ionic evaluations in each branch reached the electronic SCF criterion,
the two printed `PSTRESS` values are 457.0 kbar, and their VASP `enthalpy`
records contain the 45.7-GPa `PV` term. From the first to last evaluated step,
the two cell enthalpies fell by 25.56 and 6.08 meV, respectively (12.78 and
3.04 meV/GaN). The final 2.4-Å Ga--N coordination screens are 5/5 for the
negative seed and 4/4 for the positive seed; the corresponding maximum
atomic forces are still 2.11 and 0.57 eV/Å. Both `CONTCAR` geometries equal
the final evaluated geometries. These are **not** converged endpoints or
certified two-basin links, and a ten-step cap does not establish a strict TS.
Follow-up must start from the last evaluated geometry and keep the same
electronic and 45.7-GPa pressure contracts. VASP's [PSTRESS](https://vasp.at/wiki/PSTRESS)
and [ISIF](https://vasp.at/wiki/ISIF) specifications establish the pressure
unit and allowed degrees of freedom; do not add a second `PV` to VASP's
printed relaxation enthalpy. This diagnostic does not remove the existing
600-eV energy–stress derivative discrepancy.

### Measured narrow-tube outcome (2026-09-28)

The first nine-anchor 600-eV frozen tube passed all 18 raw static audits but
had a maximum leave-one-anchor-out error of `0.6684 meV/GaN`. A geometry-screened
second batch, Slurm array `27792651`, added nine path anchors (`1,6,8,9,16,18,19,21,22`)
and 18 more **unchanged-600-eV** VASP statics. All 18 tasks completed and their
OUTCARs passed input-hash, SCF, geometry, energy, force, and stress checks.
The combined maximum error fell to `0.3993 meV/GaN`, still above the
predeclared `0.10 meV/GaN` contour gate. Hence the narrow full-path tube is
measured but **not cleared for a smooth publication contour**.

The failure is localized, not an excuse to change `ENCUT`: the largest held-out
errors remain at images 20 (`0.399 meV/GaN`) and 5 (`0.344 meV/GaN`). The
same-contract force/stress-derived transverse slope changes sign at 5→6
(`+0.0191` to `−0.0240 eV/Å`) and 19→20 (`+0.0545` to `−0.0285 eV/Å`). At image
18 it is `+0.0886 eV/Å`, allowed by the ordinary `0.10 eV/Å` chain criterion;
a local quadratic fit to the signed statics puts its transverse minimum at
`q_perp≈−0.030 Å`, outside the sampled `±0.015 Å` strip. These correlations
do not by themselves prove a branch switch or invalidate the reported barrier.
They do show why standard barrier convergence cannot be conflated with a
sub-meV smooth landscape. The next surface step must test a wider or
conditionally relaxed coordinate chart, and/or targeted path refinement,
under the *same* 600-eV electronic inputs. It should not blindly add more
anchors to this same narrow strip or silently relax the interpolation gate.

Audited machine-readable records: `benchmarks/numerical_integrity/gan_600eV_path_tube_refinement_screen_20260928.json`,
`gan_600eV_path_tube_refinement_20260928.json`, and
`gan_600eV_path_tube_residual_analysis_20260928.json` in that directory.

### Wider-coordinate canaries and chart decision (2026-09-28)

Widening the *joint atom–strain* transverse amplitude to `±0.20 Å` was tested
at four representative anchors under exactly the same 600-eV inputs (Slurm
`27792700`). Six of eight statics completed. Both signs at image 5 failed
immediately with VASP's direct/reciprocal Bravais conflict: the direct cell
was classified base-centered monoclinic and the reciprocal simple monoclinic
for one sign, with the classes reversed for the other. The earlier empirical
single-projection geometry window had passed both, so it is demonstrably
insufficient as a general VASP preflight. Neither `SYMPREC` nor `ENCUT` was
changed to rescue these points; the joint wide tube is **not** a complete
surface. At the successful images 15, 18 and 20, the measured excess enthalpy
at `±0.20 Å` spans roughly `22–52 meV/GaN`, a physically clearer scale than
the narrow strip, but these six points cannot support an uninterrupted plot.

An alternative chart keeps the original, already-computed cell of each VCNEB
image and uses a path-orthogonal *atomic-only* transverse vector. The reaction
coordinate `s` still contains the full variable-cell transformation; the
transverse coordinate must explicitly be labeled atomic-only, not a joint
atom–strain mode. Its `±0.20 Å` canary at images 5 and 18 (Slurm `27792726`)
passed 4/4 raw VASP audits at 600 eV, including the two image-5 cases that
failed in the joint chart. However, the image-5 excess enthalpies are
`288/347 meV/GaN` and atomic forces reach `4.03 eV/Å`: that width is too
large for a useful local-mode contour there. The atomic normal is smooth
through the central barrier segment (adjacent overlap `>0.95` for images
5–22) but not over the entire chain (minimum overlap `0.147` near the final
endpoint). The defensible next design is a *piecewise/adaptive-width*
atomic-transverse chart, with smaller widths near image 5 and explicit seams
where mode identity rotates. It needs intermediate-width holdouts and a
predeclared interpolation/error gate before any publication contour.

The raw-audited canary records are
`benchmarks/numerical_integrity/gan_600eV_wide_tube_canary_20260928.json` and
`gan_600eV_atomic_tube_canary_20260928.json`; the latter derives from the
geometry-only atomic-chart definition in
`gan_600eV_atomic_tube_feasibility_20260928.json`. These studies explain a
coordinate-design problem, not a reason to change the agreed electronic
parameters or the default `0.10 eV/Å` NEB threshold.

The completed *bounded* central-segment test kept the fixed-seed atomic-only
normal over images 5–22 (minimum adjacent overlap `0.954` there) and sampled
`q_atom=±0.05 Å` at preselected anchors `5,8,11,14,17,18,20,22`, while reusing
the audited `q=0` chain. The predeclared gate for a smooth **central-segment**
off-path excess-enthalpy interpolation is a maximum leave-one-anchor-out
error of `1.0 meV/GaN`. This is about a 5% scale for a 20-meV transverse
feature and far below the roughly 0.39-eV/GaN path barrier; it is an
interpolation target, not a new DFT convergence or default NEB threshold.
The source generator and Slurm array `27792746` retain the original 600-eV
INCAR/KPOINTS/POTCAR hashes. Its gate failure calls for inspecting the
localized error and coordinate branch before deciding the next step. This central chart cannot be extended
through the endpoints without a separate mode-identity seam analysis.

**2026-09-28 central atomic-tube outcome.** Array `27792746` completed all
16 signed off-path statics; the raw audit found a maximum linear
leave-one-anchor-out (LOO) excess-enthalpy error of `6.598 meV/GaN`, so the
predeclared `1.0 meV/GaN` gate failed. A bounded, same-600-eV refinement
`27792774` added signed `q=±0.05 Å` anchors at images `6,7,19,21` and
independent `q=±0.025 Å` half-steps at images `8,20`; all 12 statics passed
raw audit. The combined 12-anchor linear LOO maximum fell to
`2.329 meV/GaN`, still a gate failure. The even transverse-curvature
half-step discrepancy was only `0.0541%` at worst, below its separately
predeclared 5% gate. Thus the local `q` curvature is resolved at those two
anchors, but the along-path interpolation is not.

An *offline, post-hoc* comparison on exactly the same audited points gives
maximum LOO errors of `2.329` (linear), `2.725` (PCHIP), `2.689` (Akima),
and `2.433 meV/GaN` (natural cubic). None meets `1.0 meV/GaN`; model
selection on these already-seen holdouts would not independently validate a
new contour even if one had passed. No cutoff, k mesh, PAW, or SCF setting
was changed. The result is a raw sampled central transverse cut, **not** a
validated smooth GaN 2D enthalpy contour. Sources and reproducible comparison
are `gan_600eV_atomic_tube_central_20260928.json`,
`gan_600eV_atomic_tube_refinement_20260928.json`, and
`gan_600eV_atomic_tube_interpolants_20260928.json` under
`benchmarks/numerical_integrity/`.

**Why interpolation fails (same data, no new DFT).** Decompose each signed
pair into an even `q²` response and an odd response at `q=0.05 Å`. The
maximum linear LOO errors are `1.907 meV/GaN` for the even part (worst at
image 8) and `1.229 meV/GaN` for the odd part (worst at image 20). Replacing
the odd interpolation with the already archived `q=0` projected force at
each held-out image still leaves a maximum signed error of
`1.935 meV/GaN`; it is a post-hoc diagnostic, not an independent pass. Near
images 19–21 the transverse force-derived slope changes from about
`-0.084` to `-0.090` to `+0.010 eV/Å`, while adjacent fixed-seed normal
overlaps are approximately `0.955` around image 20. The even excess also
has localized structure near image 8. This identifies a changing local
path/coordinate response; it does not by itself prove a discontinuous
physical branch or a software error. Analyze those image geometries and
mode transport before considering a predeclared, bounded *along-s* holdout.
Do not respond by changing `ENCUT`, relaxing the original NEB threshold, or
promoting an after-the-fact fit. The reproducible component report is
`benchmarks/numerical_integrity/gan_600eV_atomic_tube_components_20260928.json`.

An explicitly **exploratory, discrete-sample** two-panel figure now accompanies
this audit at `paper/VARNEB_CPC/figures/gan_600eV_atomic_tube_samples.*`.
Panel (a) shows the unchanged 29-image 600-eV path (`0.33849 eV/GaN`
forward barrier); panel (b) places the 18 central `q=0` images and all 28
audited off-path statics at their actual `(s,q)` coordinates. There is no
color interpolation or contour. The tracked PDF/SVG/PNG, source CSV and hash/QA
JSON accompany a reproducible local 600-dpi TIFF export (repository-ignored);
the QA record states the failed `2.329 > 1.0 meV/GaN` LOO
gate. It is suitable as a transparent interim or supplementary diagnostic,
not a substitute for the requested validated smooth whole-path 2D surface.

**Dense central-grid test submitted 2026-09-28.** To turn the barrier-region
chart into an explicitly sampled two-coordinate cut, Slurm array `27793051`
on hf `hfacnormal01` fills the central images 5–22 at
`q_atom = −0.05, −0.025, 0, +0.025, +0.05 Å`. Of the 90 chart coordinates,
18 path centers and 28 signed off-path statics are already audited and
reused; only 44 missing statics are submitted (32 MPI each, at most four
concurrent). The transverse displacement is **atomic-only** and retains each
audited variable-cell image cell. All 44 inputs have the original 600-eV
INCAR/KPOINTS/POTCAR hashes and pass the geometry preflight. Images 10 and
16 were designated as prospective along-`s` holdouts before DFT; the maximum
error of a linear-in-`s`, quadratic-in-`q` prediction of their off-path
*excess* enthalpies must be
≤1.0 meV/GaN. The inner `q=±0.025 Å` values separately test a quadratic
interpolant fitted to `q=0,±0.05 Å`, with a ≤1.0 meV/GaN analysis gate
fixed before inspecting the new outputs (unlike the along-`s` gate, this
second gate was not in the submitted manifest). The raw audit must establish
SCF termination, unchanged geometry/input hashes, finite energy/force/stress
and both validation errors
before a smooth central-segment contour is promoted. Even a pass would not
turn this atomic-only central chart into a global two-phonon PES, close the
endpoint mode-identity seams, or certify the full-variable-cell TS.

Figure contract for the next GaN composite: the hero panel plots the
path-adapted central `(s,q_atom)` enthalpy cut and every raw grid point,
with the audited `q=0` VCNEB path traced explicitly. A subordinate panel
shows the complete 29-image one-dimensional barrier profile and identifies
the central chart domain. The figure uses Python/matplotlib at double-column
width (~183 mm), editable PDF/SVG plus 600-dpi TIFF, source CSV and hash/QA
JSON. Panel labels are normal-weight `(a)`/`(b)`, all four spines are shown,
ticks and legends remain readable at final size, and any legend has a
semi-transparent white framed background away from the data. If either
validation gate fails, show measured samples and report the failure rather
than styling an unvalidated interpolant as a publication PES.

**Dense-grid outcome (same day).** All 44 array elements of `27793051`
finished `COMPLETED 0:0`; the independent raw audit passed every input hash,
SCF-completion marker, unchanged atomic/cell geometry, and finite energy,
force, and stress. The eight prospective along-`s` holdout values have a
maximum interpolation error of **0.16578 meV/GaN**, below the pre-submission
1.0 meV/GaN gate. The 36 inner-`q` tests have a maximum quadratic-model
error of **0.15888 meV/GaN**, below the separately fixed pre-audit 1.0
meV/GaN gate. The sampled transverse excess enthalpy spans 0–20.60
meV/GaN. The validated central contour and full 29-image path comparison are
`paper/VARNEB_CPC/figures/gan_600eV_atomic_dense_surface.pdf` (also SVG,
PNG and a local 600-dpi TIFF); the adjacent `_source_data.csv` has all 90
measured coordinates with source hashes, and `_qa.json` records the figure
scope. The machine-readable audit is
`benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json`.
The figure is a **validated central frozen atomic-transverse enthalpy cut**
around a variable-cell path—not a global two-phonon surface, conditional
minimum-energy surface, finite-temperature free energy, or certified TS.

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

## CPC figure contract (Python/matplotlib)

Core conclusion: the 600-eV GaN path has a local coupled atom–strain negative
direction, while a path-adapted transverse coordinate reveals how the
enthalpy changes away from the full B4→B1 VCNEB centerline. Use a full-width
quantitative composite (183 mm) with two aligned, untitled panels: (a) the
whole-path `H(s,q_perp)` tube as the hero, preserving every computed point
and drawing the actual `q_perp=0` chain; (b) the near-saddle frozen
`H(q_u,q_v)` cut with its grid and axial holdouts. Mark panels `(a)` and `(b)`
in normal-weight large type, keep all four axis spines and legible labels,
and use framed semitransparent white legends away from data. Export editable
PDF/SVG plus 600-dpi TIFF, a source-data CSV and a hash/QA record.

Only interpolate within the sampled coordinate rectangle and show raw markers.
For panel (a), interpolate the *off-path excess enthalpy* along arc length,
not the barrier itself; use the predeclared `1.0 meV/GaN` maximum LOO gate,
which the current central atomic chart does **not** pass. Do not publish a
smooth contour as validated until an independently assessed chart passes.
For panel (b), use a declared
local quadratic interpolant only if the independent half-step holdout error
passes `0.20 meV/GaN`. Neither contour is a free-energy landscape or proof of
a globally unique MEP. The transverse normal is transported and is not a
single global Γ phonon; the small q-range and numerical stress/energy
consistency limit belong in the caption, not hidden in source files.

For a stationary enthalpy minimum or an index-one saddle, the derivative
`dH/ds` along a smooth MEP is zero at that point. An index-one saddle has
negative reaction-direction curvature; stable endpoints have positive
curvature. The GaN barrier is therefore not two independent barrierless
transitions joined at a stable intermediate: the shared TS is unstable.
This statement concerns the 0 K enthalpy path, not a finite-temperature
free-energy profile or the slope of a discrete image-number plot.
