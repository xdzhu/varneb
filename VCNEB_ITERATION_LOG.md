# VC-NEB Iteration Log

## HfO2 conservative-step retry and preflight diagnostic (`fab47e8`)

- A deliberately malformed basis-directory setting was tested first. The
  calculator preflight rejected image 0 before ABACUS launched and reported the
  image, isolated work directory, launcher context, and missing basis path;
  this directory is retained as a negative validation record.
- Before the corrected retry, `cu26` was checked at `2026-09-12T06:54:53+08:00`:
  40 cores, load averages `0.04/0.20/0.20`, and no active DFT process.
- A new branch resumed the complete FIRE-40 trajectory with
  `FIRE(maxstep=0.05 A)`. The seven complete steps had maximum generalized
  forces `0.617654`, `0.622142`, `0.633117`, `0.652366`, `0.678677`,
  `0.705180`, and `0.725727 eV/A`. The force continued to increase, so the
  run was stopped early and all processes belonging to this work directory
  were cleaned up.
- This isolates two facts: the calculator error path is actionable, but a
  smaller optimizer step alone does not make the current HfO2 path convergent.
  The trial is diagnostic data only and is not a transition barrier.

Compact record: `outputs/abacus_hfo2_levelB_manifest.json`.

## HfO2 static force/stress diagnostic (`9d79e72`)

- Before the run, `cu17` was checked at `2026-09-12T07:07:47+08:00`:
  40 cores, load averages `0.08/0.03/1.97`, and no active DFT process. The
  run completed with no remaining DFT process.
- The latest complete seven-image FIRE retry trajectory was evaluated with
  zero optimizer steps at Dojo-NC-FR, 60 Ry, 1x1x1 k points, and `scf_thr=1e-6`.
  The highest image was image 3 with relative enthalpy `0.305708 eV` and
  NEB residual `0.621982 eV/A`.
- Image 2, rather than the highest-energy image, carried the largest interior
  raw atomic force (`0.674825 eV/A`) and stress component
  (`0.016611 eV/A^3`). This spatial mismatch is a concrete indicator that
  the current linear/restarted path needs a mechanism-aware path improvement
  before a production VC-NEB barrier is attempted.
- The generated summary contains cell metrics and true/spring/NEB force
  decompositions for every image. It is diagnostic data only.

Compact records: `outputs/abacus_hfo2_levelB_static_diag_cu17_manifest.json` and
`outputs/abacus_hfo2_levelB_static_diag_cu17_summary.json`.

## Release-and-refine workflow (`6f88a34`)

- Before the cluster run, `cu17` was checked at `2026-09-12T07:21:09+08:00`:
  40 cores, load averages `0.04/0.56/1.73`, and no active DFT process. The
  run used only the analytic coupled calculator and left no DFT process.
- A projected collective-mode stage converged in 12 FIRE steps with barrier
  `0.286000 eV` and maximum generalized force `0.009891 eV/A`.
- Releasing the constraint and refining the same chain in the full extended
  space converged in 37 FIRE steps with barrier `0.250044 eV` and maximum
  generalized force `0.008655 eV/A`. The barrier dropped by `0.035956 eV`.
- This validates the intended semantics: a constrained mechanism path can be
  useful for searching, but its barrier must be distinguished from the final
  unconstrained VC-NEB result.

Compact records: `outputs/release_and_refine_model/manifest.json`,
`outputs/release_and_refine_model/release_and_refine_cluster.json`, and
`outputs/release_and_refine_model/cluster_run.log`.

## CPC manuscript draft (`6f88a34`)

- Added `paper/vcneb_CPC.tex` using the requested five-part structure:
  `Introduction -> Theory -> Software -> Examples -> Conclusions/Availability`.
  No standalone Benchmarks section was introduced.
- The draft reports only supported results: analytic finite differences, CI and
  fixed-cell ASE reduction, release-and-refine behavior, calculator smoke tests,
  and the explicitly non-converged HfO2 diagnostic case.
- Local TeX Live compilation succeeded and produced a five-page PDF. The
  rendered pages were visually inspected; only minor underfull/overfull box
  warnings remain for later copy editing.

## HfO2 Level B optimizer trial (`04bf038`)

- Checked `cu17` before both launches; it had no active DFT/MPI process and 40
  cores.  The run used independently relaxed 12-atom HfO2 T/PO endpoints,
  Dojo-NC-FR, 60 Ry, 1x1x1 k points, and `scf_thr=1e-6`.
- A fresh seven-image FIRE branch completed all 40 requested steps and wrote
  41 complete snapshots.  Its final generalized force was `0.500065 eV/A`
  and its provisional enthalpy barrier was `0.394209 eV`; this is a complete
  but non-converged restart source.
- A separate LBFGS branch resumed the last complete FIRE trajectory.  Its
  force rose from `0.617654` to `1.542518 eV/A` in the first four steps, so
  the task was stopped after six complete snapshots.  No other user's process
  was touched; the remaining logs and snapshots are retained.
- The result demonstrates that optimizer choice and SCF-backed path geometry
  need explicit convergence evidence.  Neither provisional barrier is a
  publishable HfO2 transition barrier.
- Compact record: `outputs/abacus_hfo2_levelB_manifest.json`.

## Installable package skeleton (`0.1.0.dev0`)

- Added `pyproject.toml`, package version metadata, `vcneb --version`, and
  `python -m vcneb --version` entry points.
- The core dependency set is only ASE and NumPy; plotting and development
  packages are optional extras.  VASP/ABACUS remain runtime calculator
  configurations rather than mandatory package dependencies.
- Local packaging verification built
  `vcneb-0.1.0.dev0-py3-none-any.whl` successfully and both version entry
  points returned `0.1.0.dev0`.

## Fixed-cell ASE comparison (`818a87f`)

- Added `examples/run_fixed_cell_ase_comparison.py`, which sends the same
  perturbed seven-image path and analytic calculator to ASE CINEB and to
  VCNEB with a zero cell mask.
- The cluster run on `cu26` used ASE `3.23.1b1`, FIRE, and `fmax=0.01 eV/A`.
  The two paths had the same endpoint reaction energy `0.1774 eV`; barriers
  were `0.28117077 eV` (ASE) and `0.28116144 eV` (VCNEB), an absolute
  difference of `9.33e-6 eV`.
- The result confirms the fixed-cell reduction at the energy/path level for
  this analytic calculator.  It does not replace a physical DFT comparison.
- Compact record: `outputs/fixed_cell_ase_comparison.json`.

## Direction-basis conflict diagnostics (`91380bb`)

- Added `direction_basis_conflicts()` for atomic/cell direction bases
  combined with component masks.  It reports partially clipped columns, fully
  inactive columns, rank loss, and removed-component norms without hiding the
  conflict inside a later SVD failure.
- Added regression coverage for partial clipping and complete deactivation.
  The full suite passed on `cu26` after the change.

## CI saddle diagnostics (`7c58250`)

- Added `VCNEB.highest_image_index()` and `VCNEB.saddle_diagnostics()`.
  The latter reports the highest interior image, residual generalized force,
  tangent/normal residuals, and a nonuniform-spacing local tangent-curvature
  estimate.  The documentation explicitly limits this curvature to a path
  diagnostic rather than a full Hessian.
- Extended `tests/check_vcneb_forces.py` with a seven-image coupled atom/cell
  double-well test.  The cluster run on `cu26` passed the complete regression
  suite, including `climbing_image_saddle_regression=ok`.
- The analytic CI path recovered the known `0.25 eV` barrier, the saddle at
  `(q_x, strain_xx)=(0.5, 0.125)` within tolerance, a residual generalized
  force below `0.01 eV/A`, and negative local tangent curvature.

## 2026-09-12 05:00 +08:00

Focus: compare an independent fresh CI-VCNEB branch against the resumed
no-climb branch.

- Checked `cu24` immediately before launch; it had no active ABACUS/VASP/MPI
  process. Started from the relaxed-endpoint linear path in a new directory,
  with climbing enabled, so the previous no-climb run was not overwritten.
- The 7-image, 20-step ABACUS CI smoke completed with barrier `0.784444 eV`,
  reaction enthalpy `-0.747857 eV`, and final generalized force
  `0.504943 eV/A`. Its validator report was `ok` and it produced 21 complete
  snapshots.
- The no-climb branch after 30 total resumed steps had barrier `0 eV` because
  all internal images were below the higher-energy T endpoint, while its
  generalized force was still `0.390162 eV/A`. This is a non-converged path,
  not evidence of a zero physical barrier.
- The fresh CI branch retained an internal high-energy image, demonstrating
  why a climbing image and independent branch comparison are required.

Both branches remain explicitly marked as low-precision, non-converged
engineering data in `outputs/abacus_hfo2_relaxed_endpoint_manifest.json`.

## 2026-09-12 04:03 +08:00

Focus: relax the HfO2 endpoints and run a variable-cell band from the relaxed
structures.

- Checked `cu17` before endpoint work and `cu24` before the band; both had no
  active DFT task at launch. Each node has 40 cores.
- Used ASE `FrechetCellFilter` with ABACUS Dojo-NC-FR at `60 Ry`, `1x1x1`
  k points and 40 MPI ranks. Thirty-step endpoint trials reached
  `0.0383 eV/A` for T and `0.0788 eV/A` for PO at the loose `0.1 eV/A`
  target; their volumes were `143.833` and `138.589 A^3`.
- Used the two endpoint `CONTCAR` files for a 7-image, 10-step, no-climb
  VC-NEB smoke on `cu24`. Energy decreased overall, all seven image
  directories and SCF logs were produced, and the ABACUS validator passed.
- The final generalized force was `0.5130 eV/A`; the returned barrier
  `0.786010 eV` is not converged and is not a physical claim.
- The ABACUS driver now writes summary JSON/TXT and final POSCAR files, so the
  same `vcneb.traj` can be used for a later resume test.

Evidence is recorded in `outputs/abacus_hfo2_relaxed_endpoint_manifest.json`.

## 2026-09-12 03:44 +08:00

Focus: complete the corrected analytic convergence matrix.

- Checked `cu17` immediately before use; it had no active DFT process.
- The first matrix run was rejected because the convergence driver had omitted
  the toy potential's final `exx=0.25` endpoint strain, which produced a false
  `delta=0.1774 eV`.
- After fixing that endpoint, reran 54 combinations: 5/7/9 images, weak/
  medium/strong springs, `cell_scale=4/5/6 A`, and FIRE/LBFGS.
- All 54 cases converged at `fmax=0.005 eV/A`, recovered barrier `0.25 eV`,
  and recovered reaction energy `0 eV`. Full data are in
  `outputs/vcneb_convergence_matrix_fixed.json`; run metadata are in
  `outputs/vcneb_convergence_manifest.json`.

The result validates the model implementation and exposes optimizer/metric
sensitivity through force and iteration columns. It does not replace the
pending DFT image/cell/electronic convergence matrix.

## 2026-09-12 03:31 +08:00

Focus: validate the physical ABACUS path before production VC-NEB.

Execution:

- Checked `cu17` immediately before each calculation; it had 40 cores, no
  active VASP/ABACUS job and low load at the single-image and multi-image
  launch points.
- Ran the HfO2 12-atom T endpoint with ABACUS Dojo-NC-FR, 40 MPI ranks,
  `ecutwfc=60 Ry`, `1x1x1` k points, and `cal_force=cal_stress=out_stru=1`.
- The single-image calculation returned energy, finite stress and atomic
  forces; the generated inputs passed the static ABACUS validator.
- Ran a seven-image, two-step VC-NEB smoke from the mapped T endpoint to the
  mapped PO endpoint. Every image produced its own ABACUS work directory and
  the driver completed with an output band and trajectory.

Evidence:

- Results and node metadata are stored in `outputs/abacus_hfo2_smoke_manifest.json`.
- Single-image smoke passed. The multi-image run is an engineering smoke only:
  final optimizer `fmax` was about `1.805086 eV/A`, so the barrier
  `0.677095 eV` is explicitly not a converged physical result.
- The ABACUS driver now writes `vcneb_summary.json`, `vcneb_summary.txt`,
  final POSCAR files and a band plot for future runs.

## 2026-09-12 03:28 +08:00

Focus: run the updated analytic model cases on the shared cluster.

Execution:

- Checked `cu17` immediately before use: 40 cores, load average
  `0.08, 0.06, 0.05`, no visible user jobs.
- Ran `examples/run_toy_vcneb.py` and
  `examples/run_hfo2_t_po_model_vcneb.py` from the shared project directory
  using the ICU Python environment.
- Toy model result: barrier `0.250004 eV`, endpoint difference `0.000000 eV`.
- HfO2 12-atom T->PO synthetic model result: barrier `0.800023 eV`, endpoint
  difference approximately zero.
- Saved the result and limitation statement in
  `outputs/vcneb_model_case_manifest.json`.

Verification:

- The same node then passed the complete regression suite, including the
  non-diagonal/pressure cell-force checks, mode constraints, mask validation,
  and cell-validity guards; `compileall` also passed.

Scientific boundary:

- These are analytic/synthetic model cases. They verify end-to-end mechanics,
  not a physical VASP/ABACUS phase-transition barrier. Production DFT remains
  a separate P3/P5 task.

## 2026-09-12 03:20 +08:00

Focus: strengthen generalized cell-force validation and multi-mode handling.

Changes:

- Added a rotation-invariant metric-cell analytic calculator to the regression
  suite, so non-diagonal deformation and non-zero pressure are tested against
  a physically symmetric stress tensor.
- Added pressure/enthalpy finite-difference diagnostics that separate the zero-
  pressure stress contribution from the `P dV` contribution.
- Added least-squares coefficients for `project_path_onto_modes()` when the
  supplied modes are not orthogonal; array and iterable mode inputs are both
  accepted.
- Added finite 0/1 validation for `atom_mask` and `cell_mask`.

Verification:

- On the shared `cu17` environment, after a fresh occupancy check, the full
  regression suite passed with
  `max_zero_pressure_cell_force_error=3.156e-11`,
  `max_pressure_volume_force_error=8.380e-10`, and
  `max_pressure_cell_force_error=8.343e-10`.
- The mode-subspace, projected-mode, non-orthogonal-mode and mask-validation
  regressions also passed.

Interpretation:

- The earlier `1.546e-02` discrepancy was a test-model defect: the toy model
  supplied a non-symmetric stress for an energy depending on rotational cell
  components, while ASE exposes the physical symmetric stress tensor. The
  core cell-force formula was not changed; the corrected physical test now
  validates it to numerical precision.

## 2026-09-12 03:09 +08:00

Focus: implement and document explicit mode-constrained VCNEB semantics.

Changes:

- Added `build_mode_basis()` for atomic/cell mode vectors in the extended
  VCNEB coordinate space.
- Added `build_direction_basis()` for per-atom directions and cell-direction
  basis vectors.
- Added `mode_basis` and `constraint_mode` to `VCNEB` and `run_vcneb`.
  `subspace` validates endpoint compatibility and projects interior images and
  updates into one affine subspace; `projected` preserves each initial image's
  non-mode component while projecting its updates.
- Added regression checks for strict subspace endpoint rejection, force/update
  projection and direction-basis layout.
- Updated the README, theory draft and project plan to distinguish guided
  initial paths, projected dynamics and strict subspace paths.

Verification:

- Local Python syntax and `git diff --check` passed.
- After checking `cu17` occupancy (load average `0.00, 0.01, 0.05`, 40 cores,
  no visible user jobs), synchronized commit `36493a0` to the shared project
  directory and ran the full regression suite.
- All existing checks plus `mode_subspace_constraint_regression=ok` and
  `projected_mode_constraint_regression=ok` passed.

Remaining work:

- Verify strict mode projection against a reduced-coordinate implementation and
  finite differences.
- Add sparse projector/conflict diagnostics and release-then-refine workflow.

## 2026-09-12 02:58 +08:00

Manuscript organization decision:

- The VCNEB paper will use the compact structure
  `Introduction -> Theory -> Software -> Examples -> Conclusions/Availability`.
- Verification, convergence studies and literature/reference-implementation
  comparisons will be presented as subsections, tables or figures within
  `Examples`, rather than as a standalone `Benchmarks` section.

## 2026-09-12 02:45 +08:00

Focus: establish the GitHub-synced project baseline and validate it on the
shared cluster filesystem.

Changes:

- Initialized the local research directory as a Git worktree and connected it
  to `https://github.com/xdzhu/vcneb`.
- Pushed baseline commit `31b1b40` to `origin/main`.
- Added `.gitignore` rules that keep DFT restart/output files, caches and
  generated workspaces out of the source repository.
- Added `VCNEB_PROJECT_PLAN.md` with P0--P7 milestones and acceptance criteria.
- Added `outputs/vcneb_p0_baseline_manifest.json` with the node, environment,
  commit and regression results.
- Confirmed the CPC reference package at
  `D:\Work\Zstar\zstar-article\submission_packages\ZStar_CPC_pdflatex` and its
  `elsarticle` structure for the future VCNEB manuscript.

Cluster verification:

- Checked `cu17` before use: 40 cores, load average approximately
  `0.00, 0.01, 0.05`, no user jobs visible.
- Synced the tracked baseline to the shared path
  `/home/zhuxd/abacus/agent-runs/20260912-vcneb-p0`.
- Ran `/home/zhuxd/Software/anaconda3/envs/icu/bin/python3
  tests/check_vcneb_forces.py` on `cu17`; all 11 regression checks passed.
- The node image does not provide `git`, so source synchronization used a
  `git archive` generated from commit `31b1b40`. The project nodes share the
  same filesystem and program environment; future runs need only one shared
  source sync plus a fresh pre-run occupancy check.

Status:

- P0 baseline freeze: source and regression baseline established.
- P1 theory specification: next critical task, starting with cell-force/stress
  finite-difference conventions and the `cell_scale` metric.

## 2026-09-12

Focus: USPEX 10.6 implementation forensics and mode-guided VC-NEB support.

Findings:

- `H:\\ReSearch\\VCNEB\\USPEX\\USPEX_v10.6.tar.gz` is a MATLAB Runtime self-extracting ELF distribution, not a readable source archive.
- The public Qian paper and USPEX manual remain the clean algorithm references; no USPEX binary code was copied or reverse engineered.
- ABINIT provides `iatfix/iatfixx/iatfixy/iatfixz` component constraints and `nconeq/iatcon/wtatcon` linear force constraints, but not a turnkey arbitrary phonon-mode NEB.

Changes:

- Added `vcneb/modes.py` with `Mode`, `mode_guided_path`, and `project_path_onto_modes`.
- Added `atom_mask` to `VCNEB` and `run_vcneb` for ABINIT-like component constraints.
- Added regression coverage for endpoint-preserving mode seeds, modal projection, and inactive atomic components.
- Added `outputs/uspex_vcneb_mode_analysis.md` and updated `README_VCNEB.md` with the OpenVCNEB positioning and mode semantics.

Important semantic boundary:

- `mode_guided_path` is an initial-path generator; subsequent unconstrained NEB still targets the full-space MEP.
- `atom_mask` creates a constrained calculation and must be reported as such.
- Strict arbitrary mode-subspace MEPs remain a separate feature requiring endpoint-subspace validation.

Verification:

- Local `py_compile` static check passed for the modified Python modules.
- Full regression execution is reserved for `235 -> cu05` per the cluster-only calculation policy.

## 2026-06-21 08:29 +08:00

Focus: harden ABACUS static validation against disabled required VC-NEB outputs.

Changes:

- Updated `scripts/validate_vcneb_inputs.py`.
  - Added ABACUS flag parsing for numeric and common boolean forms.
  - `cal_force`, `cal_stress`, and `out_stru` must now be enabled, not just present.
  - Templates that set any of those keys to `0`, `false`, `.false.`, `off`, or similar disabled values now fail the dry static check before a VC-NEB run is launched.
- Extended `tests/check_vcneb_forces.py`.
  - Renamed and broadened the ABACUS validator regression to `abacus_validator_required_flags_regression=ok`.
  - The test still rejects missing `out_stru`, accepts a minimal valid template, and now rejects disabled `cal_force`, `cal_stress`, and `out_stru` one by one.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `abacus_validator_required_flags_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.

Remote checks:

- Synced lightweight source/tests into `/home/zhuxd/abacus/agent-runs/20260621-vcneb` on `235`.
- On `235`:
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_validator_required_flags_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`:
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_validator_required_flags_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The next hardening target remains a tiny one-step real ABACUS force/stress smoke in isolated scratch, using copied input assets only.
- Resume still restores geometry only, not optimizer internal state.

## 2026-06-19 06:40 +08:00

Focus: first automation run for variable-cell NEB hardening.

Changes:

- Hardened `cell_mask` handling in `vcneb/core.py`.
  - Inactive cell components are now preserved from each current image during optimizer `set_x` updates instead of being reset to the initial endpoint cell.
  - NEB tangent, spring, climbing-image, reaction-coordinate, and final generalized forces now exclude masked cell degrees of freedom.
- Added validation for `dynamic_relaxation` and `dynamic_energy_scale` constructor inputs.
- Extended `tests/check_vcneb_forces.py` with a regression check that masked cell components neither move through `set_x` nor receive spring/true-force leakage.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.
- Lightweight adapter checks:
  - `vcneb.vasp.read_vasp_input_params("initial_state/relax")`: pass, parsed 19 INCAR keys and KPOINTS keys `gamma`, `kpts`, `reciprocal`.
  - `vcneb.abacus.make_ase_abacus_factory(...)`: pass for factory construction. Local ASE does not include `ase.calculators.abacus`.

Remote checks:

- Copied only `README_VCNEB.md`, `vcneb/`, `examples/`, and `tests/` to `/home/zhuxd/abacus/agent-runs/20260619-vcneb` on `235`.
- On `235` with conda env `icu`, ASE `3.23.1b1`:
  - `python tests/check_vcneb_forces.py`: pass.
  - `python examples/run_toy_vcneb.py`: pass.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`, same scratch directory and ASE `3.23.1b1`:
  - `python tests/check_vcneb_forces.py`: pass.
  - `python examples/run_toy_vcneb.py`: pass.
  - `python -m compileall vcneb examples tests`: pass.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- Local ASE is `3.28.0`; remote `235/cu05` ASE is `3.23.1b1`. The optimizer compatibility smoke passes on both, but the ABACUS calculator adapter still needs validation in an actual minimal ABACUS force/stress call before production use.
- Next hardening target: add restart/resume support for VC-NEB trajectories and a dry-run validator for VASP/ABACUS per-image work directories.

## 2026-06-19 10:55 +08:00

Focus: real 12-atom crystal endpoint fixture for HfO2 tetragonal-to-polar-orthorhombic VC-NEB validation.

Changes:

- Added `scripts/setup_hfo2_t_po_validation.py`.
  - Builds a 12-atom `sqrt(2) x sqrt(2) x 1` tetragonal HfO2 cell from the local VASP T-phase endpoint.
  - Parses the local ABACUS O/PO-phase `STRU`, maps atoms species-wise with Hungarian assignment and a translation search, and writes a portable fixture under `validation/hfo2_t_to_po/`.
  - Avoids ASE-version-dependent `make_supercell` ordering by enumerating the tetragonal supercell explicitly.
- Added `examples/run_hfo2_t_po_model_vcneb.py`.
  - Uses the real mapped HfO2 atom/cell path with a synthetic endpoint double-well calculator.
  - This is a geometry, mapping, optimizer, MIC, and cell-DOF smoke test, not a physical HfO2 potential.
- Added `scripts/validate_vcneb_inputs.py`.
  - Performs static dry-run checks for VASP template inputs plus per-image POSCAR-like files.
  - Performs static dry-run checks for ABACUS `INPUT/KPT/STRU`, force/stress keys, and referenced pseudo/orbital files when present.
  - Does not launch VASP or ABACUS.
- Fixed a `mic=True` interpolation bug in `vcneb/core.py`.
  - The interpolation was correctly choosing the shortest fractional endpoint, but then reset the final image back to the wrapped input endpoint.
  - This reintroduced a periodic discontinuity and produced a spurious huge barrier on the HfO2 path.
- Inspected the earlier remote implementation under `/home/zhuxd/abacus/8.dielec/vcneb`.
  - Reused the useful design direction: fractional coordinates plus cell degrees of freedom, ASE optimizer compatibility, species-wise Hungarian/MIC atom mapping, and VASP per-image directories.
  - Kept the newer implementation's stricter choices: endpoints are excluded from optimizer DOFs, cell force conversion follows the ASE UnitCellFilter row-vector convention, and force/cell transforms have finite-difference regression tests.

Endpoint validation:

- `validation/hfo2_t_to_po/T_HfO2_12.vasp`: `P4_2/nmc (137)` by `spglib`.
- `validation/hfo2_t_to_po/PO_HfO2_12_mapped.vasp`: `Pca2_1 (29)` by `spglib`.
- Mapping translation: approximately `[0.0, 0.5, 0.5]`.
- Atom mapping initial index to original final index: `[1, 0, 2, 3, 4, 7, 6, 5, 8, 11, 10, 9]`.
- Mean assignment distance: `0.5176093344 A`.
- Initial path segment distances: `0.361440, 0.360895, 0.360359, 0.359832, 0.359315, 0.358807 A`.

Local checks:

- `python scripts/setup_hfo2_t_po_validation.py`: pass.
- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
- `python examples/run_hfo2_t_po_model_vcneb.py`: pass.
  - `barrier_eV=0.800023`
  - `delta_eV=-0.000000`
  - `enthalpies_eV=0.000000 0.261003 0.633792 0.800023 0.620812 0.248664 0.000000`
- `python scripts/validate_vcneb_inputs.py --mode vasp --template-dir initial_state/relax --image-root validation/hfo2_t_to_po`: pass.
  - `status=ok`
  - `n_images=7`
  - `atom_counts=12 12 12 12 12 12 12`
- `python -m compileall vcneb examples scripts tests`: pass.

Remote checks:

- Synced `README_VCNEB.md`, `VCNEB_ITERATION_LOG.md`, `vcneb/`, `examples/`, `scripts/`, `tests/`, and `validation/` to `/home/zhuxd/abacus/agent-runs/20260618-vcneb` on `235`.
- On `cu05` through `ssh 235`:
  - `python3 scripts/setup_hfo2_t_po_validation.py`: pass, same mapping and path distances as local.
  - `python3 examples/run_hfo2_t_po_model_vcneb.py`: pass.
    - `barrier_eV=0.800023`
    - `delta_eV=0.000000`
    - `enthalpies_eV=0.000000 0.261003 0.633792 0.800023 0.620812 0.248664 0.000000`
  - `python3 tests/check_vcneb_forces.py`: pass.
    - `max_fractional_force_error=3.178e-11`
    - `max_cell_force_error=5.657e-11`
    - `cell_mask_regression=ok`
  - `python3 scripts/validate_vcneb_inputs.py --mode vasp --template-dir /home/zhuxd/abacus/8.dielec/vcneb/vasp_images/img_00 --image-root validation/hfo2_t_to_po`: pass.
    - `status=ok`
    - `n_images=7`
    - `atom_counts=12 12 12 12 12 12 12`
  - `python3 -m compileall vcneb examples scripts tests`: pass.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The HfO2 model test confirms the VC-NEB coordinate/force machinery on a real mapped crystal path, but the next validation should run at least a very cheap one-step force/stress call for every image with either VASP or ABACUS.
- Next hardening target: generate a tiny per-image DFT smoke case from the validated HfO2 path and run one force/stress call per image before any full VC-NEB optimization.

## 2026-06-19 12:36 +08:00

Focus: restart/resume support for interrupted calculator-backed VC-NEB runs.

Changes:

- Added `read_chain_trajectory(...)` in `vcneb/core.py`.
  - Reads the latest complete `n_images` frame group from the flat ASE trajectory written by `run_vcneb`.
  - Ignores an interrupted partial tail instead of treating it as a valid final chain.
  - Validates atom order and 3D cells before returning the chain.
- Added `apply_chain_state(...)` in `vcneb/core.py`.
  - Copies resumed cells and positions into existing image objects while preserving attached calculators.
- Added `trajectory_mode` to `run_vcneb(...)` with supported modes `w` and `a`.
- Wired `--resume` and `--resume-trajectory` into both `examples/run_vcneb_vasp.py` and `examples/run_vcneb_abacus.py`.
  - When resuming from the active `workdir/vcneb.traj`, new optimizer frames append to the existing trajectory.
- Documented the resume contract in `README_VCNEB.md`.
- Extended `tests/check_vcneb_forces.py` with a trajectory-resume regression covering complete snapshots, partial-tail recovery, and calculator preservation.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.
- Local append-mode smoke using `run_vcneb(..., trajectory_mode="a")`: pass.
  - `append_resume_smoke=ok 4`

Remote checks:

- Copied only `README_VCNEB.md`, `vcneb/`, `examples/`, and `tests/` to `/home/zhuxd/abacus/agent-runs/20260619-vcneb` on `235`.
- On `235` with the lightweight Python/ASE environment:
  - `python tests/check_vcneb_forces.py`: pass, including `trajectory_resume_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`, same scratch directory:
  - `python tests/check_vcneb_forces.py`: pass, including `trajectory_resume_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The shell printed `未找到 PID ... 进程的工作目录信息。` before the Python output; it did not affect the smoke checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- Resume restores chain geometry only; optimizer internal state is not restored, so resumed runs restart the optimizer dynamics from the resumed coordinates.
- Next hardening target remains a tiny one-step force/stress DFT smoke per image, preferably ABACUS first because its ASE adapter/version surface is less stable than VASP's.

## 2026-06-19 18:40 +08:00

Focus: harden append-mode resume snapshots so interrupted calculator-backed runs do not lose earlier per-step snapshots.

Changes:

- Updated `run_vcneb(...)` in `vcneb/core.py`.
  - Added `snapshot_start` for explicit caller control.
  - In `trajectory_mode="a"`, the default snapshot counter now continues from the next available `step_####` or `chain_step_####.traj` index under `snapshot_dir`.
  - Existing snapshot files are no longer overwritten during append-mode resume.
- Added a focused regression to `tests/check_vcneb_forces.py`.
  - Seeds `snapshots/chain_step_0000.traj` with a marker.
  - Runs a one-step append-mode VC-NEB smoke.
  - Verifies the marker survives and `chain_step_0001.traj` is created.
- Documented append-mode snapshot numbering in `README_VCNEB.md`.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.

Remote checks:

- Synced only lightweight changed files to `/home/zhuxd/abacus/agent-runs/20260619-vcneb` on `235`.
- On `235`:
  - `python tests/check_vcneb_forces.py`: pass, including `snapshot_append_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - `ase.calculators.abacus.Abacus` imports with ASE `3.23.1b1`, signature `(profile=None, directory='.', **kwargs)`.
- On `cu05` through `ssh 235`:
  - `python tests/check_vcneb_forces.py`: pass, including `snapshot_append_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- Resume still restores chain geometry only; optimizer internal state is not restored.
- Next hardening target remains a tiny one-step force/stress DFT smoke per image, preferably ABACUS first, using an isolated scratch case and copied input assets only.

## 2026-06-20 02:14 +08:00

Focus: harden the ASE ABACUS calculator factory against version-specific constructor behavior.

Changes:

- Updated `vcneb/abacus.py`.
  - `make_ase_abacus_factory(...)` now rejects ambiguous `command` plus `profile` inputs.
  - When `command` is provided, it is converted to `ase.calculators.abacus.AbacusProfile(command)` and passed as `profile=...`.
  - This avoids leaking `command` into ABACUS native input parameters on ASE builds whose `Abacus` signature is `(profile=None, directory='.', **kwargs)`.
  - Required VC-NEB force/stress keys are still enforced: `cal_force=1`, `cal_stress=1`, and `out_stru=1`.
- Extended `tests/check_vcneb_forces.py`.
  - Added a mocked `ase.calculators.abacus` regression that verifies command-to-profile conversion, required key defaults, extra kwargs preservation, and rejection of simultaneous `command`/`profile`.

Local checks:

- Local Python/ASE probe:
  - Python `3.10.9`.
  - ASE `3.28.0`.
  - Local ASE does not include `ase.calculators.abacus`, so the ABACUS regression uses a fake module.
- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
  - `abacus_command_profile_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.

Remote checks:

- Copied only `README_VCNEB.md`, `vcneb/`, `examples/`, and `tests/` to `/home/zhuxd/abacus/agent-runs/20260620-vcneb` on `235`.
- On `235`:
  - ASE `3.23.1b1`.
  - `ase.calculators.abacus.Abacus` imports with signature `(profile=None, directory='.', **kwargs)`.
  - `AbacusProfile` imports with signature `(command, pseudo_dir=None, basis_dir=None, **kwargs)`.
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_command_profile_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`:
  - ASE `3.23.1b1`.
  - `ase.calculators.abacus.Abacus` imports with signature `(profile=None, directory='.', **kwargs)`.
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_command_profile_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output on `cu05` and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The ABACUS adapter now matches the discovered ASE constructor/profile API on `235/cu05`, but still needs a tiny one-step real ABACUS force/stress smoke in scratch before production VC-NEB use.
- Resume still restores chain geometry only; optimizer internal state is not restored.

## 2026-06-20 08:22 +08:00

Focus: harden ASE optimizer API shape compatibility for non-FIRE optimizers.

Changes:

- Updated `vcneb/core.py`.
  - `VCNEB.get_masses()` now returns one mass per optimizer coordinate row (`len(chain)`), matching `get_positions().shape[0]`.
  - The previous value returned one mass per scalar DOF, e.g. 24 masses for 8 optimizer rows in a 4-image one-atom band.
- Extended `tests/check_vcneb_forces.py`.
  - Added `optimizer_api_shapes_regression=ok`.
  - The regression verifies `__len__`, `get_positions()`, and `get_masses()` agree.
  - It also runs one lightweight analytic VC-NEB step through both `BFGS` and `LBFGS`.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
  - `optimizer_api_shapes_regression=ok`
  - `abacus_command_profile_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.
- `python examples/run_hfo2_t_po_model_vcneb.py`: pass.
  - `barrier_eV=0.800023`
  - `delta_eV=-0.000000`

Remote checks:

- Copied only lightweight source/docs/tests into `/home/zhuxd/abacus/agent-runs/20260620-vcneb` on `235`.
- On `235`:
  - Python `3.10.9`, ASE `3.23.1b1`.
  - `python tests/check_vcneb_forces.py`: pass, including `optimizer_api_shapes_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`:
  - `python tests/check_vcneb_forces.py`: pass, including `optimizer_api_shapes_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The next hardening target remains a tiny one-step real ABACUS force/stress smoke in isolated scratch, because the adapter constructor is now covered but the full external calculator call is still not.
- Resume still restores chain geometry only; optimizer internal state is not restored.

## 2026-06-20 14:21 +08:00

Focus: harden cached energy/force state after ASE optimizer coordinate updates.

Changes:

- Updated `vcneb/core.py`.
  - `VCNEB.set_x(...)` now invalidates cached enthalpies and generalized forces after successfully applying optimizer coordinates.
  - This prevents `get_value()`, `enthalpies`, and `barrier()` from reporting stale values when queried after `set_positions(...)` / `set_x(...)` but before a fresh force evaluation.
- Extended `tests/check_vcneb_forces.py`.
  - Added `set_x_cache_invalidation_regression=ok`.
  - The regression populates the enthalpy cache, moves the active image through `set_x(...)`, and verifies the next energy/barrier query reflects the new image geometry.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
  - `optimizer_api_shapes_regression=ok`
  - `set_x_cache_invalidation_regression=ok`
  - `abacus_command_profile_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.
- `python examples/run_hfo2_t_po_model_vcneb.py`: pass.
  - `barrier_eV=0.800023`
  - `delta_eV=-0.000000`

Remote checks:

- Copied only lightweight source/docs/tests into `/home/zhuxd/abacus/agent-runs/20260620-vcneb` on `235`.
- On `235`:
  - Python `3.10.9`, ASE `3.23.1b1`.
  - `python tests/check_vcneb_forces.py`: pass, including `set_x_cache_invalidation_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`:
  - Python `3.10.9`, ASE `3.23.1b1`.
  - `python tests/check_vcneb_forces.py`: pass, including `set_x_cache_invalidation_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The next hardening target remains a tiny one-step real ABACUS force/stress smoke in isolated scratch, using copied input assets only.
- Resume still restores chain geometry only; optimizer internal state is not restored.

## 2026-06-21 02:32 +08:00

Focus: harden ABACUS adapter and static validation so required VC-NEB outputs cannot be accidentally disabled.

Changes:

- Updated `vcneb/abacus.py`.
  - `make_ase_abacus_factory(...)` now merges user parameters and calculator kwargs first, then force-enforces `cal_force=1`, `cal_stress=1`, and `out_stru=1`.
  - This closes the previous gap where `setdefault(...)` allowed accidental `0` values to pass through, leaving VC-NEB without force, stress, or output-structure data.
- Updated `scripts/validate_vcneb_inputs.py`.
  - ABACUS template validation now checks for `out_stru` alongside `basis_type`, `cal_force`, and `cal_stress`.
- Extended `tests/check_vcneb_forces.py`.
  - The mocked ABACUS factory regression now verifies required keys override conflicting user inputs.
  - Added `abacus_validator_out_stru_regression=ok`, which rejects a minimal ABACUS template missing `out_stru` and accepts it after adding the key.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
  - `optimizer_api_shapes_regression=ok`
  - `set_x_cache_invalidation_regression=ok`
  - `parallel_endpoint_ownership_regression=ok`
  - `abacus_command_profile_regression=ok`
  - `abacus_validator_out_stru_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.

Remote checks:

- Synced lightweight source/docs/tests into `/home/zhuxd/abacus/agent-runs/20260621-vcneb` on `235`.
- On `235`:
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_validator_out_stru_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - Real ASE ABACUS constructor-only probe: pass; conflicting `cal_force`, `cal_stress`, and `out_stru` inputs all landed as `1` on the constructed `Abacus` calculator with `AbacusProfile`.
- On `cu05` through `ssh 235`:
  - `python tests/check_vcneb_forces.py`: pass, including `abacus_validator_out_stru_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The next hardening target remains a tiny one-step real ABACUS force/stress smoke in isolated scratch, using copied input assets only.
- Resume still restores geometry only, not optimizer internal state.

## 2026-06-20 20:23 +08:00

Focus: harden parallel image ownership so MPI-style VC-NEB force reductions do not double-count endpoint enthalpies.

Changes:

- Updated `vcneb/core.py`.
  - Endpoint images are now assigned to rank 0 in the ownership table.
  - In `parallel=True` mode, only rank 0 evaluates endpoint images before the collective sum.
  - Interior image ownership remains round-robin over active images.
  - This prevents endpoint enthalpies from being counted once per rank when `_enthalpy_and_force(...)` reduces values through ASE's parallel world.
- Extended `tests/check_vcneb_forces.py`.
  - Added `parallel_endpoint_ownership_regression=ok`.
  - The regression monkeypatches a lightweight four-rank `world` object and verifies endpoint ownership is rank-0-only while interiors keep the expected round-robin assignment.

Local checks:

- `python tests/check_vcneb_forces.py`: pass.
  - `max_fractional_force_error=3.178e-11`
  - `max_cell_force_error=5.657e-11`
  - `cell_mask_regression=ok`
  - `trajectory_resume_regression=ok`
  - `snapshot_append_regression=ok`
  - `optimizer_api_shapes_regression=ok`
  - `set_x_cache_invalidation_regression=ok`
  - `parallel_endpoint_ownership_regression=ok`
  - `abacus_command_profile_regression=ok`
- `python examples/run_toy_vcneb.py`: pass.
  - `barrier_eV=0.250004`
  - `delta_eV=0.000000`
- `python -m compileall vcneb examples tests`: pass.
- `python examples/run_hfo2_t_po_model_vcneb.py`: pass.
  - `barrier_eV=0.800023`
  - `delta_eV=-0.000000`

Remote checks:

- Synced only lightweight source/docs/tests into `/home/zhuxd/abacus/agent-runs/20260620-vcneb` on `235`.
- On `235`:
  - `python tests/check_vcneb_forces.py`: pass, including `parallel_endpoint_ownership_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
- On `cu05` through `ssh 235`:
  - `python tests/check_vcneb_forces.py`: pass, including `parallel_endpoint_ownership_regression=ok`.
  - `python examples/run_toy_vcneb.py`: pass, `barrier_eV=0.250004`, `delta_eV=0.000000`.
  - `python -m compileall vcneb examples tests`: pass.
  - The known Chinese process-directory warning printed before Python output and did not affect the checks.

Remaining limitations / next actions:

- No production VASP or ABACUS DFT was launched in this run.
- The next hardening target remains a tiny one-step real ABACUS force/stress smoke in isolated scratch, using copied input assets only.
- Resume still restores chain geometry only; optimizer internal state is not restored.
# 2026-09-12

Focus: USPEX 10.6 implementation forensics and mode-guided VC-NEB support.

## VASP single-image smoke (05:28, `b470f52`)

- Before the run, `cu26` was checked through gateway `235`: 40 cores, load average
  `0.07/0.04/0.05`, and no active ABACUS/VASP/MPI process.
- Synced the repository at commit `b470f52` into the shared run directory and used
  the legacy GaN `POSCAR/POTCAR` from `/home/zhuxd/abacus/8.dielec/vcneb/vasp`.
- Generated a minimal static VASP template with `IBRION=-1`, `NSW=0`, `ISIF=2`,
  `ISYM=0`, `EDIFF=1E-5`, and a Gamma-centered `2x2x2` mesh.  The run used
  `mpirun -np 40 .../vasp_std` with VASP 6.3.2.
- ASE calculator preflight reported energy, forces, stress, variable-cell support,
  an image-local directory, and a unique directory; all issues were empty.
- Result: 4 atoms, energy `-22.57971864 eV`, maximum force
  `0.19039522 eV/A`, finite stress matrix, volume `45.72832351 A^3`.
- `scripts/validate_vcneb_inputs.py --mode vasp` returned `status=ok` for the
  template and generated image directory.  The result is an interface smoke only,
  not a converged GaN or HfO2 VC-NEB barrier.
- Compact record: `outputs/vasp_gan_single_image_manifest.json`.

Findings:

- `H:\\ReSearch\\VCNEB\\USPEX\\USPEX_v10.6.tar.gz` is a MATLAB Runtime self-extracting ELF distribution, not a readable source archive.
- The public Qian paper and USPEX manual remain the clean algorithm references; no USPEX binary code was copied or reverse engineered.
- ABINIT provides `iatfix/iatfixx/iatfixy/iatfixz` component constraints and `nconeq/iatcon/wtatcon` linear force constraints, but not a turnkey arbitrary phonon-mode NEB.

Changes:

- Added `vcneb/modes.py` with `Mode`, `mode_guided_path`, and `project_path_onto_modes`.
- Added `atom_mask` to `VCNEB` and `run_vcneb` for ABINIT-like component constraints.
- Added regression coverage for endpoint-preserving mode seeds, modal projection, and inactive atomic components.
- Added `outputs/uspex_vcneb_mode_analysis.md` and updated `README_VCNEB.md` with the OpenVCNEB positioning and mode semantics.

Important semantic boundary:

- `mode_guided_path` is an initial-path generator; subsequent unconstrained NEB still targets the full-space MEP.
- `atom_mask` creates a constrained calculation and must be reported as such.
- Strict arbitrary mode-subspace MEPs remain a separate feature requiring endpoint-subspace validation.

Verification:

- Local `py_compile` static check passed for the modified Python modules.
- Full regression execution is reserved for `235 -> cu05` per the cluster-only calculation policy.
