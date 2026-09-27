# VARNEB manuscript evidence gate — 2026-09-27

This is a claim gate, not a request to change a calculator parameter or to run
CI. The current draft PDF has 12 pages; the intended CPC article is about 10.
The working tree contains unrelated and ongoing user research files: stage only
the specific artifacts being revised.

## What can already be shown

1. **BTO T→C:** seven real ABACUS VCNEB images and a cubic Γ-mode projection
   are in `figures/bto_gamma_mode_path.*`. The main-text exploratory figure
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
3. **GaN B4→B1, 45.7 GPa:** `figures/gan_multibackend_validation.*` shows
   five converged *one-dimensional* variable-cell enthalpy paths, normalized
   per GaN: ABACUS/VASP/QE/ABINIT/CP2K forward barriers
   0.3274/0.3385/0.3297/0.2924/0.2928 eV/GaN. This is reasonably near
   Qian et al.'s ~0.34 eV/GaN tetragonal-route reference while calculator
   protocols are not identical. `figures/gan_joint_mode_coupling.*` shows a
   local atom–strain instability candidate; its 1000-eV VASP Hessian and
   600-eV endpoint Γ optical bases are only compared geometrically. Neither
   a unique full-variable-cell TS nor a GaN two-dimensional PES is certified.

## Minimum submission gates for a focused CPC software paper

- Keep the 59-point BTO contour explicitly labeled as an exploratory frozen
  cut; leave the conditional surface as staged work until its stronger gate
  below is complete. The main text now follows this narrower claim.
- Freeze figure-to-data provenance and captions: explicit `E` versus
  `H=E+PV`, pressure, formula-unit normalization, reference zero, number of
  computed points, interpolation status, and backend-specific settings.
- Reconcile code, examples, quick-start README, detailed manual, and paper
  around the same backend/optimizer interface and supported feature matrix.
  Re-run relevant unit/integration checks locally; do not trigger CI merely
  by pushing ordinary commits.
- Tighten the 12-page draft toward about 10 pages without deleting method
  derivation or evidence limits; check PDF layout, bibliography, figure
  legibility, source-data availability, and journal formatting.

## Stronger mode-surface / TS claims: evidence still needed

- **BTO conditional 2D surface:** the second curvature step and direct
  soft-direction checks at both `Q_z=0.9` sites are complete. Three
  independently audited static starts at `(0.75,0.15)` and their
  conditional relaxations (`27787487`) are complete. An independent audit of
  all 77 DFT points finds unique gradient/stress-eligible cache points in the
  `Q_y=0,+,−` basins; the pre-declared low-energy center prediction misses by
  only `0.476 meV/BTO`, under its `2 meV/BTO` local line. The frozen runner
  did not explicitly serialize all branch terminal records; a cache-only
  deterministic replay reconstructed all three with zero new DFT calls.
  Local holdout curvature remains unchecked. Competing
  branch and continuity checks across the sampled domain; quantified
  interpolation error and a clear definition of released atomic/strain
  variables. Confirm that the actual T endpoint is represented by, or is
  explicitly off, the chosen two-mode manifold. Only then promote a contour
  from sampled candidates to a conditional surface.
- **GaN true variable-cell TS mode:** refine the highest image to a
  stationary point under one energy protocol at 45.7 GPa; verify the complete
  atom-plus-cell Hessian has exactly one unstable direction, appropriate
  endpoint connections, pressure/metric conventions, and finite-difference
  convergence. The current negative local joint direction is evidence for
  coupling, not full TS certification. A 2D GaN landscape is a separate
  computational project, not implied by its converged 1D path.
- **Follow-on rather than CPC prerequisites:** bilayer hBN fixed-cell sliding
  versus controlled in-plane cell relaxation, and HfO₂ T→PO multi-mode
  analysis. These require distinct boundary/functional/dispersion and mode
  contracts; they should enter a higher-level methods paper only after a
  clean reproducible result exists.
