# VC-NEB Iteration Log

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
