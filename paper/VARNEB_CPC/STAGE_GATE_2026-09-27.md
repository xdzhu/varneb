# VARNEB manuscript evidence gate — 2026-09-27

This is a claim gate, not a request to change a calculator parameter or to run
CI. The 2026-09-28 two-column `elsarticle` draft now compiles to 10 main-text
pages plus a separate two-page supplement. The BTO restricted-sheet figure,
GaN five-backend path comparison, GaN local atom--strain diagnostic, and
validated central GaN two-coordinate cut remain in the main text. The
BTO/HfO2 material-control figure and all-image GaN endpoint-mode profile are
Supplementary Figs. S1 and S2; their original source data remain linked.
This meets the length target without shrinking the plotted labels. The working
deadline is an October 3 data freeze and October 10
submission-ready manuscript package; do not make HfO₂/2D extensions or a strict
GaN TS certificate prerequisites for this focused CPC paper.
The working tree contains unrelated and ongoing user research files: stage only
the specific artifacts being revised.

## What can already be shown

1. **BTO T→C:** seven real ABACUS VCNEB images and a cubic Γ-mode projection
   are in `figures/bto_gamma_mode_path.*`. The archival exploratory figure
   `figures/bto_frozen_soft_mode_landscape.*`, copied from the independently
   audited 59-point source under `outputs/batio3_t_to_c_pbe100_dzp10au/`,
   contains 59 audited *frozen cubic-cell* DFT samples. Its smooth contour is
   display interpolation (35 interior leave-one-out points; maximum absolute
   error 5.14 meV/BTO on the final set; earlier independent center and edge
   holdouts reached 32.57 and 14.83 meV/BTO). The variable-cell path is
   projected onto that plane
   but leaves it: the contour is neither a relaxed conditional PES nor a
   T→C barrier. The cubic Γ source is the five-atom `1×1×1` cell, with
   ABACUS PBE/`ecutwfc=100 Ry`/10 au DZP/`4×4×4` electronic k points.
2. **BTO conditional mode plane, stage figure:**
   `figures/bto_conditional_four_point_stage_2026-09-27.*` shows four
   independently calculated fixed-`(Q_z,Q_x)` coordinates and the lower
   `±Q_y` branches. No interpolation is drawn. Selected energies relative to
   the cubic reference are −102.20, −104.97, −108.55, and −109.93 meV/BTO
   at `(0.6,0)`, `(0.6,0.3)`, `(0.9,0)`, `(0.9,0.3)` in
   `sqrt(amu) Å`. The branch is 57.47, 49.98, 34.46, and 28.31 meV/BTO below
   the separately relaxed `Q_y=0` branch. At `Q_z=0.9`, the raw-audited
   16D orthogonal Hessians at two finite-difference steps have positive lowest
   eigenvalues: 0.005699/0.005817 at 0.05 and 0.005720/0.005844 at 0.10,
   in eV/(amu Å²). The respective lowest directions overlap above 0.999999.
   Four additional signed DFT probes per site give positive direct
   energy- and force-derived curvatures along those directions. All 145/146
   cached DFT points, respectively, passed an independent raw-output audit.
   This strengthens the *local numerical screen*, but does not certify a
   global minimum, branch continuity, or a continuous conditional PES.
   A separate preselected center `(0.75,0.15)` is now audited with the same
   ABACUS/PBE/100-Ry/10-au-DZP contract. After mirror-branch alignment, its
   energy differs from the four-corner bilinear prediction by 0.476 meV/BTO,
   while the atomic-plus-strain coordinate discrepancy is
   `0.1254 sqrt(amu) Å`, predominantly strain. All five selected low-energy
   branches pass the declared local gradient/stress and two-step 16D
   curvature screens, but the weakest curvature lacks a rigorous numerical
   error bound. This five-point pilot does not promote the four-point stage
   figure to a continuous conditional-PES contour. Source:
   `benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json`.
   Separately, a *symmetry-restricted* `Q_y=0` variable-cell sheet now uses
   nine raw-audited training nodes, a sixth-order even-mode model, and two
   genuine prospective DFT holdouts under the same ABACUS 100-Ry/10-au-DZP
   contract. Both full raw-output audits passed (13 and 18 DFT evaluations).
   The independent energy errors are −0.150 and −0.977 meV/BTO against a
   preregistered ±2-meV/BTO gate; atom–strain seed, orthogonal-gradient, and
   2-kbar stress gates also pass. The T endpoint projects to
   `(Q_z,Q_x,Q_y)=(1.2043,0,0) sqrt(amu) Å`, with nonzero stable-mode
   displacement and strain, and C is the origin. The new main-text
   `figures/bto_qy0_restricted_even_mode_sheet.*` is drawn only inside the
   measured-coordinate hull and marks all DFT points. It is *not* an
   unrestricted globally minimized PES or finite-temperature FES. Source:
   `benchmarks/numerical_integrity/bto_qy0_even_mode_two_holdout_gate_20260928.json`.
3. **GaN B4→B1, 45.7 GPa:** `figures/gan_multibackend_validation.*` shows
   five converged *one-dimensional* variable-cell enthalpy paths, normalized
   per GaN: ABACUS/VASP/QE/ABINIT/CP2K forward barriers
   0.3274/0.3385/0.3297/0.2924/0.2928 eV/GaN. This is reasonably near
   Qian et al.'s ~0.34 eV/GaN tetragonal-route reference while calculator
   protocols are not identical. `figures/gan_joint_mode_600eV.*` shows a
   local atom–strain instability candidate under the **same original 600-eV
   VASP contract** as the path and endpoint Γ bases. Its `0.0232 eV/Å`
   energy–force mismatch still prevents a strict full-variable-cell TS
   certificate. The separate local `q_u/q_v` frozen quadratic cut now has
   eight same-input grid statics and four axial half-step holdouts, with
   maximum holdout error `0.04974 < 0.20 meV/GaN`; its contour is confined
   to those measured coordinate bounds. This does not validate off-axis
   interior points, an orthogonally relaxed PES, or the whole path. The
   supplementary candidate
   `figures/gan_600eV_atomic_tube_samples.*` plots 28 audited central
   off-path statics without interpolation; the `2.329 meV/GaN` maximum LOO
   error exceeds the `1.0` gate. A subsequent **dense central-segment**
   atomic-transverse cut at images 5–22 has 90 measured coordinates, reusing
   18 path centers and 28 earlier signed statics and adding 44 raw-audited
   points. The prospective along-path holdout error is
   `0.16578 < 1.0 meV/GaN`; the inner-transverse quadratic check is
   `0.15888 < 1.0 meV/GaN`. Its validated central contour is
   `figures/gan_600eV_atomic_dense_surface.*` with source CSV and QA JSON.
   This remains a frozen atomic-only transverse chart around a variable-cell
   centerline, not a global two-phonon plane, all-path surface, conditional
   minimum, or TS certificate. Earlier 1000-eV diagnostics are archival and
   excluded.

## Minimum submission gates for a focused CPC software paper

- Keep the 59-point frozen-cell BTO contour archival. The main text now uses
  the blinded-check-passing `Q_y=0` restricted model contour with all measured
  locations visible; never call it an unrestricted conditional PES. No further
  grid densification is required merely for a smoother illustration.
- Freeze figure-to-data provenance and captions: explicit `E` versus
  `H=E+PV`, pressure, formula-unit normalization, reference zero, number of
  computed points, interpolation status, and backend-specific settings.
- Reconcile code, examples, quick-start README, detailed manual, and paper
  around the same backend/optimizer interface and supported feature matrix.
  Re-run relevant unit/integration checks locally; do not trigger CI merely
  by pushing ordinary commits. The 2026-09-28 tracked local suite has
  549 passes and one skip, documented in
  `docs/VARNEB_LOCAL_REGRESSION_2026-09-28.md`; it does not replace remote
  material-output audits or a clean release-build check.
- Preserve the now 10-page main text and two-page supplement while checking
  journal formatting, final PDF layout, bibliography, figure legibility, and
  source-data availability. Do not regain length by dropping method derivation
  or evidence limits. Author order,
  affiliations, CRediT roles, funding and the competing-interest declaration
  still contain explicit draft placeholders and require author confirmation.

## Stronger mode-surface / TS claims: evidence still needed

- **BTO unrestricted conditional 2D surface (not a CPC prerequisite):** the second curvature step and direct
  soft-direction checks at both `Q_z=0.9` sites are complete. Three
  independently audited static starts at `(0.75,0.15)` and their
  conditional relaxations (`27787487`) are complete. An independent audit of
  all 77 DFT points finds unique gradient/stress-eligible cache points in the
  `Q_y=0,+,−` basins; the pre-declared low-energy center prediction misses by
  only `0.476 meV/BTO`, under its `2 meV/BTO` local line. The frozen runner
  did not explicitly serialize all branch terminal records; a cache-only
  deterministic replay reconstructed all three with zero new DFT calls.
  Local holdout curvature has since been screened at two finite-difference
  steps, but its weak positive value is not a rigorously certified minimum.
  The new restricted `Q_y=0` model has a clear released-coordinate definition
  and two prospective pointwise errors under 1 meV/BTO, but competing branch,
  continuity, and restricted-curvature checks across its entire domain remain
  open. The actual T endpoint is on the restricted coordinate slice, not on
  the frozen two-eigenvector plane. Do not upgrade the local model contour to
  a globally certified conditional PES on pointwise holdouts alone.
  There is also an endpoint-level obstruction to the proposed *fully
  minimized* two-coordinate lower envelope: cubic C has three independent
  unstable Γ directions, but fixed `(Q_z,Q_x)` leaves `Q_y` free. At C, the
  free `Q_y` direction has negative curvature, so the lower envelope cannot
  contain C at its reference energy. The source-bound proof is
  `benchmarks/numerical_integrity/bto_soft_triplet_conditional_endpoint_20260928.json`.
  Do not launch a dense lower-envelope grid as a surrogate T→C surface;
  explicitly choose a symmetry-restricted `Q_y=0` sheet, a three-soft-mode
  surface, or a path-adapted chart before further sampling.
- **GaN true variable-cell TS mode:** refine the highest image to a
  stationary point under one energy protocol at 45.7 GPa; verify the complete
  atom-plus-cell Hessian has exactly one unstable direction, appropriate
  endpoint connections, pressure/metric conventions, and finite-difference
  convergence. The current negative local joint direction is evidence for
  coupling, not full TS certification. A validated **central** two-coordinate
  cut now exists, but endpoint seams and a full-path chart remain separate
  work. The two-sided native VASP 600-eV/45.7-GPa basin-relaxation array
  `27793292` completed and its full raw audit verifies identical first ten
  signed-pilot evaluations, complete SCFs and B4/B1 structural returns.
  With the user-accepted 2-kbar stress gate, native B4 passes force/stress,
  while native B1 has 2.528 kbar residual and narrowly fails stress. Fresh
  statics on both final geometries reproduce the archived endpoint enthalpies
  within 0.207/0.026 meV per four-atom cell but have larger 8.165/3.726
  kbar raw-stress residuals. Native/static `PSTRESS` is 457/0 kbar; VASP's
  pressure-shifted `external pressure` line is not the common stress gate.
  These data support two-basin identification, **not**
  a stress-certified two-basin link or an index-one variable-cell TS. Do not
  mix native constant-basis relaxation enthalpies with fresh-static/path
  enthalpies; keep the same 600-eV electronic contract and diagnose the
  stress/basis-history discrepancy before promoting the claim.
- **Follow-on rather than CPC prerequisites:** bilayer hBN fixed-cell sliding
  versus controlled in-plane cell relaxation, and HfO₂ T→PO multi-mode
  analysis. These require distinct boundary/functional/dispersion and mode
  contracts; they should enter a higher-level methods paper only after a
  clean reproducible result exists.
