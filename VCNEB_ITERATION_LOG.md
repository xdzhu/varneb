# VC-NEB Iteration Log

## VARNEB 0.0.1 release workflow (`2026-09-13`)

- `pyproject.toml` and `vcneb/version.py` now report `0.0.1`; the wheel builds
  as `varneb-0.0.1-py3-none-any.whl`.
- Added `.github/workflows/publish-pypi.yml` using GitHub OIDC Trusted
  Publishing. It runs only for an explicit `v*` tag, a published Release, or a
  manually confirmed `publish=true` dispatch; ordinary pushes do not publish.
- Commit `35c0dbb` and tag `v0.0.1` were pushed to `xdzhu/varneb`. PyPI
  visibility remains pending the repository's Trusted Publisher registration.
- The first tag-triggered run (`34763435900`) built successfully but the publish
  job returned `invalid-publisher`. The OIDC claims confirm the exact expected
  repository/workflow/environment tuple: `xdzhu/varneb`,
  `.github/workflows/publish-pypi.yml`, `pypi`. After correcting the pending
  publisher, rerun the failed publish job or dispatch this workflow manually.
- After the pending publisher was added, rerunning `34763435900` completed
  successfully. `pip index versions varneb` and a clean `pip download` now
  resolve `varneb==0.0.1` from PyPI.

## VARNEB repository rename (`2026-09-13`)

- Upstream repository is now `https://github.com/xdzhu/varneb`; local `origin`
  fetch/push URLs were updated without pushing or starting CI.
- The distribution/project name in `pyproject.toml` is now `varneb`.  The
  `varneb` console entry point is primary, while the `vcneb` console alias and
  Python import package remain for backwards compatibility.

## HfO₂ parallel parser race fixed and resumed (`2026-09-13`)

- Initial ordinary parallel VCNEB job `27677945` completed three ordinary
  steps (`fmax=0.966334 -> 0.800934 eV/A`) before failing in image 3 with an
  empty force array.  ABACUS SCF itself completed; the traceback and worker
  manifest are retained remotely under the same workdir, so this is not a
  physical path failure.
- Root cause was the installed ASE ABACUS adapter writing the process-global
  `ase_sort.dat` while four image threads were concurrently preparing inputs.
  The adapter now suppresses that file for identity species ordering (the
  Hf4O8 fixtures are grouped Hf then O), while preserving upstream behavior
  for non-identity serial calculations.  A regression check passes on the
  remote ICU Python environment.
- The stale global sort file was moved aside, and job `27678004` resumed from
  the complete step-2 chain snapshot with the same 128-task/4x32-MPI/no-CI
  settings.  Its first resumed evaluation is `fmax=0.800934 eV/A` and is
  running without the previous parser warning.
- At the latest poll the resumed job has reached ordinary steps 0--6 with
  `fmax=0.800934, 0.738536, 0.677679, 0.626060, 0.591285, 0.559697,
  0.514160 eV/A`; all worker
  batches remain `status=ok` and the Slurm stderr is empty.

## HfO₂ endpoint gate passed; ordinary VCNEB launched (`2026-09-13`)

- PO continuation jobs `27677874` and `27677918` completed normally from the
  same trajectory.  The final PO summary passes the strict gate with
  `max_generalized_force=0.0004973 eV/A` and `max_abs_stress=0.07076 kbar`;
  T `27677491` passes with `0.0004988 eV/A` and `0.01841 kbar`.  Both promoted
  structures are 12-atom `Hf4O8` conventional-cell endpoints; the promotion
  report is `validation/hfo2_t_to_po/endpoint_promotion_gate.json` on `hf`.
- A calculator-free preflight on the promoted endpoints passed for a 7-image
  log-strain path: minimum distance `2.02496 A`, maximum deformation `0.04841`,
  no folded junctions, and all seven ABACUS calculator reports expose energy,
  forces, stress, variable-cell support, and isolated directories.  The local
  copy is `outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_preflight_hfo2/`.
- Ordinary HfO₂ VCNEB job `27677945` is now RUNNING on `hfacnormal01` with one
  128-task controller allocation and four concurrent 32-MPI image workers.
  It uses 100 Ry, Hf/O Orb-DZP-10au, 2x2x2 k points, FIRE, `k=0.2`, no climb,
  and a 300-step budget.  Step 0/1/2 residuals are `0.966334`, `0.881405`,
  and `0.800934 eV/A`; this is an ordinary-path preconvergence, not a final
  barrier claim.

## HfO₂ endpoint fast-BFGS status (`2026-09-13`, jobs 27677491/27677467)

- The T endpoint continuation `27677491` completed from the preserved trajectory
  with ASE BFGS (`maxstep=0.005`, 100 Ry, Orb-DZP-10au, 32 MPI).  It passes the
  common gate: maximum generalized force `0.00049882 eV/A`, maximum stress
  `0.01841 kbar`, 12 atoms (`Hf4O8`), `converged=true`.
- The active job `27677467` is the PO endpoint, not the high-symmetry T endpoint.
  It is the lower-symmetry 12-atom orthorhombic conventional cell and is being
  relaxed in the same variable-cell space.  At the latest poll (20:10 runtime)
  it reached BFGS step 9 with monotonic `fmax=0.483720 -> 0.294759 eV/A` and
  remains RUNNING on `node148` with 32 MPI.  No rebound or numerical instability
  is present; promotion and VCNEB submission remain gated on its force and
  stress summary.
- The accepted T trajectory, summary, final structure, and Slurm metadata are
  archived under `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_relax_T_ase_bfgs_maxstep005_job27677491/`.

## BTO duplicate native-endpoint submission cancelled (`2026-09-13`)

- BTO already has the complete `100 Ry + Orb-DZP-10au` endpoint/path validation matrix: 5/7/9-image ordinary VCNEB, reverse direction, rebound-window continuation, and 6×6×6/1e-9 precision comparison.
- A native ABACUS `cell-relax` endpoint rerun (`27676929` cubic, `27676930` tetragonal) was mistakenly submitted while migrating the endpoint methodology. `27676929` failed during startup and `27676930` was cancelled; neither is production evidence and no BTO VCNEB path will be rerun.
- The native `calculation cell-relax`/BFGS migration remains required for the HfO₂ endpoint gate, where the existing external-optimizer branches did not provide a production-quality endpoint.
- PO resume `27676977` and T resume `27676976` each reached an ABACUS trust-radius breakdown after preserving `STRU_ION_D`; this is an optimizer termination, not structure divergence. The saved geometries are being continued with ABACUS `bfgs_trad` (`27677116` PO, `27677117` T).
- Endpoint optimizer choice is methodological, not a different target: with identical energy/cell/ionic degrees of freedom and thresholds, converged ASE and ABACUS-native relaxations target the same stationary minimum. The BTO FIRE endpoints are therefore valid; BFGS is the preferred solver for new HfO₂ production attempts, not a reason to invalidate or repeat BTO.
- The native-driver restart parser now handles ABACUS `STRU_ION_D` labels/magnetic flags and converts dimensionless lattice vectors using `LATTICE_CONSTANT` (Bohr) to ASE Angstrom units; regression coverage is in `tests/check_native_abacus_driver.py`.

## HfO₂ native endpoint continuation (`2026-09-13`, jobs 27677116/27677117)

- The T and PO endpoints are continuing with ABACUS-native `calculation
  cell-relax`, `relax_method=bfgs_trad`, 100 Ry, Hf/O Orb-DZP-10au, 2×2×2
  k-points, and 32 MPI per endpoint on `hfacnormal01`.
- Both Slurm jobs remain `RUNNING` with active CPU.  Native cell-relax has
  completed several cell cycles; pressure rebounds (T approximately −28→+89
  kbar, then decreasing; PO approximately −17→+30 kbar) are being observed as
  multi-step optimizer windows, not treated as physical divergence or a stop
  condition.  No endpoint has yet passed the 0.02 eV/Å and 0.1 kbar gates.
- The future ordinary HfO₂ VCNEB controller/worker template is now
  `cluster/hf_hfo2_vcneb_parallel.slurm` (default four isolated image workers
  × 32 MPI in one 128-task allocation); it is synced and syntax-checked but
  not submitted before the endpoint gate.
- Native summaries now parse both ABACUS `Largest gradient in force` records and the `TOTAL-FORCE`/`TOTAL-STRESS` blocks emitted by `bfgs_trad`, retaining the last measured maximum force and stress for the endpoint gate.

## HfO₂ native `bfgs_trad` stability audit (`2026-09-13`, jobs 27677116/27677117)

- The two native jobs were observed across multiple complete cell-relax cycles,
  rather than being stopped at the first force/pressure rebound.  They were then
  cancelled safely after a clear numerical instability developed in the PO cell.
  Full `OUT.ABACUS`, `STRU_ION_D`, stdout/stderr, and summaries are archived in
  `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_*_native_bfgs_trad_resume2_cancelled_job*/`.
- PO ended with `converged=false`, `fmax=7.55946 eV/A`,
  `max_stress=2856.90 kbar`, `pressure=1644.58 kbar`, and volume `98.05 A^3`.
  T ended with `converged=false`, `fmax=0.132648 eV/A`,
  `max_stress=27.40 kbar`, and `pressure=-26.68 kbar`.  Both retained the
  required 12-atom `Hf4O8` composition, so the failure is an optimizer/cell-step
  stability issue, not evidence for a different physical endpoint or a mismatch
  between T and PO cells.
- Endpoint method equivalence remains the governing rule: external ASE BFGS and
  ABACUS-native `cell-relax` are alternative numerical routes to the same
  stationary minimum when they use the same energy surface, cell/ionic degrees of
  freedom, constraints, and force/stress tolerances.  The next HfO₂ attempt may
  therefore use a more conservative external BFGS (or damped native restart),
  with a hard force-and-stress gate before any VCNEB submission.  BTO remains
  frozen; no duplicate BTO endpoint or path calculation is permitted.

## HfO₂ conservative external-BFGS endpoint fallback (`2026-09-13`, jobs 27677331/27677332)

- After checking `hfacnormal01` availability and an empty user queue, the
  conservative fallback `cluster/hf_hfo2_endpoint_ase.slurm` was submitted as
  jobs `27677331` (PO) and `27677332` (T), each with 32 MPI, 100 Ry,
  Orb-DZP-10au, 2×2×2 k points, ASE `BFGS`, `maxstep=0.0005`, and a 300-step
  budget. The driver now records `force_converged`, `stress_converged`, and a
  combined `converged` flag using the common `0.02 eV/A` + `0.1 kbar` gate.
- Two earlier submissions (`27677322/27677323` and `27677328/27677329`) exited
  before ABACUS because of a remote script-sync path and an O-orbital filename
  typo; they produced no DFT steps and are excluded from evidence. The corrected
  pair is running with valid ABACUS `INPUT/KPT/STRU` and copied Hf/O assets.
- The ASE endpoint driver now supports `--resume`/`--resume-step`, reopening the
  latest complete `relax.traj` frame and appending to the existing trajectory;
  this keeps endpoint continuation consistent with the VCNEB chain-resume
  protocol.

## HfO₂ maxstep revision (`2026-09-13`, jobs 27677331/27677332 → 27677466/27677467)

- The `maxstep=0.0005` branch was intentionally stopped after a verified
  multi-step window: T reached step 19 (`0.020087 eV/A`) and PO reached step 15
  (`0.483720 eV/A`) with monotonic decrease but an impractically small cell
  update. Full trajectories are archived under the corresponding
  `endpoint_relax_*maxstep0005_cancelled_job*` directories.
- The latest frames were copied without re-interpolation into jobs `27677466`
  (T) and `27677467` (PO), using ASE BFGS `maxstep=0.005`, `fmax=0.01`, 120
  steps, the same 100 Ry/Orb-DZP-10au/2×2×2/32-MPI settings. This is a step
  control change only; the endpoint target and promotion gate are unchanged.
- The reusable Slurm fallback now defaults to `maxstep=0.005`; the earlier
  `0.0005` value is retained only as an archived diagnostic, not as the
  production default.
- Added `scripts/promote_hfo2_endpoints.py` plus a regression test. It validates
  both endpoint summaries against the common force/stress gate and atomically
  publishes only passing `CONTCAR` files to the production `relaxed_T/` and
  `relaxed_PO/` directories.

## HfO₂ production-precision endpoint gate (`2026-09-13`)

- The legacy `27674272` result remains diagnostic only (60 Ry).  New endpoint
  work uses 100 Ry, Orb-DZP-10au, 2×2×2 k points, and 32 MPI ranks.
- T endpoint FIRE `27676731` showed five consecutive rising generalized-force
  values (`0.5572 → 0.6284 eV/Å`) and was cancelled only after the observation
  window; this diagnoses an unstable variable-cell optimizer setting, not a
  physical claim that the structure itself is divergent.  Its partial trajectory
  is archived locally.
- A constrained-step LBFGS retry `27676767` likewise showed five rising values
  (`0.5572 → 0.5829 eV/Å`) and was cancelled after the same observation rule;
  this is also an optimizer/cell-coupling diagnostic.  Its partial trajectory is
  archived locally.
- PO endpoint FIRE `27676732` and fixed-cell T pre-relaxation `27676802` were
  cancelled on request after preserving their partial trajectories (PO reached
  step 32 at `0.2998 eV/Å`; fixed-cell T reached step 6 at `0.2213 eV/Å`).
- Direct ASE variable-cell BFGS runs `27676827`/`27676828` were cancelled after
  preserving their diagnostics.  The endpoint workflow now uses ABACUS-native
  `calculation cell-relax` with `relax_method bfgs`: jobs `27676868` (T) and
  `27676869` (PO) generate native `OUT.ABACUS/STRU_ION_D` and a converted
  `CONTCAR`; no HfO₂ VCNEB path will be submitted until both pass the force gate.

## BaTiO3 T→C production branch reset (`2026-09-13`)

- Corrected the case-selection decision: the earlier HfO2-versus-BaTiO3 cost
  comparison was not a fair precision comparison because the HfO2 diagnostic
  branch used `60 Ry`/`1x1x1`/`1e-6` only as a low-cost Level-A/Level-B trial.
  No HfO2 production claim will be made from that setting.
- The active material branch is now BaTiO3 tetragonal→cubic.  The cluster
  templates use the same ABACUS Dojo-NC-FR setup for endpoints and images:
  `ecutwfc=100 Ry`, `4x4x4` k points, `scf_thr=1e-8`, and the explicit
  `Dojo-NC-FR/Orb-DZP-10au` directory with `Ba/Ti/O` 10 au DZP orbitals.
- The BTO VCNEB launcher now accepts `DIRECTION` (default
  `tetragonal_to_cubic`) and `ENDPOINT_TAG`, and writes isolated work
  directories so the old mixed-basis runs are not resumed or overwritten.
- Fresh 100 Ry/10 au DZP endpoint relaxations completed as jobs `27675330`
  (cubic, `node40`) and `27675331` (tetragonal, `node41`) on `hfacnormal01`.
  The ordinary seven-image T→C no-climb run is `27675388` on `node40`; CI is
  intentionally disabled until the ordinary path is stable.
- Subsequent production templates now request 32 MPI ranks (one CPU per rank)
  on `hfacnormal01`; already-running jobs retain their original allocation.
- The ABACUS launcher now writes an atomic preflight/summary record with
  calculator capabilities, per-image directory ownership, Slurm metadata and
  Git revision; `scripts/audit_vcneb_result.py` provides a calculator-free
  post-run gate.  `path_diagnostics()` now embeds final-path geometry metrics.
- Chain snapshots and per-image POSCAR snapshots now use same-directory
  temporary files plus `os.replace`, so an interrupted write cannot masquerade
  as a complete recovery point.
- The first BTO FIRE branch (`27675388`) was stopped after a force rebound;
  the complete step-12 chain was preserved.  External-trajectory recovery
  (`27675509`, 32 MPI, `maxstep=0.002`) reproduced the state without
  re-interpolation, then settled near `0.05 eV/A`.  LBFGS continuation
  (`27675598`) also plateaued/rebounded, so both branches remain diagnostic,
  not converged MEP claims.
- A 32-MPI static audit (`27675643`) of the best preserved chain reports a
  monotonic T→C enthalpy rise of `0.08712890 eV`, no interior barrier, and
  `fmax=0.057946 eV/A`; geometry remains valid (`d_min=1.8015 A`, deformation
  `0.05168`, minimum segment cosine `0.9407`).  CI is therefore correctly
  withheld pending a mechanism/path change rather than forced onto the highest
  image.
- A further 32-MPI FIRE continuation (`27675704`, `maxstep=0.003`) was started
  from the best `maxstep=0.005` snapshot.  It improved only from
  `0.029613` to `0.029365 eV/A`, then rebounded to `0.031914 eV/A`; the run was
  safely cancelled after step 2 with `chain_step_0000..0002.traj` and the full
  trajectory preserved under
  `outputs/batio3_t_to_c_pbe100_dzp10au/vcneb_t_to_c_n7_fire_maxstep003_from_step2_cancelled/`.
  This confirms a local optimizer plateau rather than a reason to enable CI.
- The first image-count matrix was then submitted with the same BTO production
  calculator and no climb: `27675759` (5 images, `node182`) and `27675760`
  (9 images, `node383`), both 32 MPI / 59168M.  At the latest checkpoint,
  5-image step 9 has `fmax=0.393455 eV/A` and 9-image step 4 has
  `fmax=0.494433 eV/A`; both remain `RUNNING` with only successful ABACUS
  image tasks.  These are active convergence branches, not yet final results.
- Added the opt-in `ThreadedCalculatorExecutor` image backend and the
  `cluster/hf_batio3_vcneb_parallel.slurm` controller template.  A real BTO
  five-image one-step smoke (`27675909`) used 128 allocated tasks and launched
  four concurrent 32-MPI exclusive job steps; all image evaluations completed
  with exit code 0 and the controller wrote a normal summary/trajectory.  The
  archived result is under
  `outputs/batio3_t_to_c_pbe100_dzp10au/vcneb_t_to_c_n5_parallel_smoke_job27675909/`.
  Its force audit is intentionally failed (`0.521268 eV/A` after one step), but
  the image-level parallel execution and barrierless profile are validated.
- The executor now supports optional per-image retry (`--image-retries`) and
  tracks attempt counts in memory; a fail-once calculator regression verifies
  that only the failed image is retried and that the successful result is retained.
- The original 9-image serial branch was safely stopped at complete step 9
  (`0.438473 eV/A`) and resumed through the parallel controller as `27675981`
  (128-task allocation, four 32-MPI workers).  It has reached parallel step 16
  at `0.259057 eV/A` with successful image steps; the original serial trajectory
  remains archived in its own work directory.
- The 5-image branch was similarly resumed from serial step 23 as `27676003`.
  Its best parallel snapshot is step 11 at `0.0388197 eV/A`, followed by a
  rebound to `0.040979 eV/A`; it was safely cancelled and statically audited
  as job `27676030`.  The audit confirms the same monotonic barrierless profile
  (`0.08712890 eV`) and valid geometry, but the force target is not met.  Both
  serial branches are copied under the local BTO output archive for recovery
  provenance.
- Operational correction: a single FIRE rebound is not a stop condition. Future
  ordinary-NEB branches use a 3--5 complete-step rebound observation window;
  only persistent rebound/plateau without a new best snapshot triggers safe
  cancellation and static audit.
- The correction was exercised immediately in BTO: `27676180` resumed the
  five-image best complete snapshot (step 11) for 12 steps with four 32-MPI
  workers.  It rebounded at step 3 (`0.036197 eV/A`) but then descended to
  `0.0231289 eV/A` at step 12, demonstrating that the earlier cancellation was
  premature.  The run remains below the `0.02` target but is archived for the
  rebound-window evidence.
- Nine-image job `27675981` completed its requested 30-step first segment at
  `0.0739323 eV/A`; it was not labeled converged and was resumed as `27676207`
  from the complete step-30 trajectory.  Both branches use the BTO 100-Ry,
  DZP-10au calculator and four concurrent 32-MPI image workers.
- The five-image continuation `27676251` then crossed the ordinary-path force
  gate at step 7 (`0.0191050 eV/A`) after the rebound-window evidence.  Its
  calculator-free audit is `status=ok`; the 5-vs-7 comparison remains
  `barrierless-consistent` with identical `0.0871289 eV` reaction enthalpy.
- Seven-image continuation `27676310` crossed the same ordinary force gate at
  step 17 (`0.0196710 eV/A`) after a real multi-step rebound window.  Its audit
  is `status=ok`, with the same `0.0871289 eV` barrierless profile and valid
  geometry.  Thus both 5- and 7-image ordinary paths now meet the force target.
- Nine-image continuation `27676207` completed another 30 steps at
  `0.0289920 eV/A` and was archived as a non-converged segment; `27676513` is
  the next continuation from its complete step-60 trajectory.  The 5/7/9
  calculator-free comparisons remain barrierless-consistent.
- Reverse C→T job `27676299` reached step 20 at `0.208863 eV/A`; its continuation
  `27676485` is being observed without early cancellation after a later rebound
  (`0.043058 -> 0.049381 eV/A`).
- Reverse continuation `27676485` subsequently crossed the force gate at step 34
  (`0.0195577 eV/A`) after the rebound window.  Audit is `status=ok`, with zero
  forward barrier, reaction enthalpy `-0.0872206 eV`, and no interior barrier.
- The complete BTO ordinary image-count matrix is now converged: five-image
  `27676251` (`0.0191050 eV/A`), seven-image `27676310` (`0.0196710 eV/A`),
  and nine-image `27676513` (`0.0199074 eV/A`).  The calculator-free comparison
  reports `barrierless-consistent`, with identical `0.0871289 eV` reaction
  enthalpy and no geometry issues.  The CI gate is therefore explicitly
  withheld and recorded in `bto_ci_gate.json` rather than forcing a climb on a
  monotonic path.
- Electronic-precision audit: `27676621` (100 Ry, 6×6×6, SCF 1e-9, zero-step
  static path) and endpoint relaxations `27676663`/`27676664` agree on a reaction
  energy of about `0.0724 eV` (path `0.0723553 eV`, endpoints `0.0724417 eV`).
  This differs from the 4×4×4 production value by about `0.01476 eV`, so the
  barrierless topology is stable but final absolute energetics should use the
  tightened electronic setting.  Record: `bto_precision_sensitivity.json`.

## HfO2 Level B continuation audit (`27674272`, 2026-09-13)

- The resumed seven-image ordinary VC-NEB branch completed all 100 requested
  FIRE steps with `climb_after=100` (CI never activated).  Slurm reports
  `COMPLETED`, exit code `0`, partition `hfacnormal01`, 8 tasks on `node109`,
  elapsed `03:12:04`; the exact submission used `FMAX=0.05`,
  `SPRING=0.2`, `MAXSTEP=0.01`, `RESUME=1`, 60 Ry, 1x1x1 k points and
  `scf_thr=1e-6`.
- The complete remote result was archived locally under
  `outputs/hfo2_vcneb_levelA_no_climb_v1_job27674272/`, including the final
  `vcneb.traj`, all retrieved chain snapshots, optimizer log, summary and
  barrier plot.  The final trajectory contains 142 complete seven-image
  chains (994 ASE frames), so this is a valid continuation rather than a
  re-interpolated path.
- The path remains non-converged: final maximum generalized force is
  `0.436268 eV/A` versus the `0.05 eV/A` target.  The provisional forward
  enthalpy barrier is `0.619562 eV`, reaction enthalpy is `-0.743338 eV`, and
  image 3 is the local interior peak.  Its residual generalized force is
  `0.215015 eV/A` and perpendicular force is `0.466786 eV/A`.
- The ordinary path is not yet stable enough for CI.  The optimizer decreases
  `fmax` to `0.312842 eV/A` at continuation step 39, then reverses and ends at
  `0.436268 eV/A`; this is consistent with the earlier mechanism/path
  mismatch, not convergence.  Final image minimum MIC distance is `1.9480 A`
  (above the `1.6 A` guard), but adjacent final extended-coordinate cosines
  are `[-0.0751, 0.9454, 0.7933, 0.8526, 0.5920]`, so the first junction is
  folded and the path is not a clean MEP.  Final cells remain nonsingular and
  the largest endpoint-referenced deformation is about `0.1277`.
- Decision: do not enable CI or submit another long DFT run from this chain.
  Next work is mechanism-aware path construction/diagnostics (and, if needed,
  a calculator-free or static-force test) before a fresh ordinary VC-NEB
  production attempt.  This result must not be reported as a physical HfO2
  transition barrier.

## ABACUS high-precision BaTiO3 staged-CI diagnostic (`27673595`)

- After the algorithm-first path preflight, ran a seven-image BaTiO3
  cubic-to-tetragonal ABACUS branch on `hf` with 100 Ry, `4x4x4` k points,
  `SCF_THR=1e-8`, `FIRE(maxstep=0.005)`, and `climb_after=15`.
- Occupancy was checked before submission. The run used 16 tasks on `node77`
  because of the high plane-wave/k-point cost; task count remains a per-run
  resource choice rather than a fixed 40-core policy.
- The ordinary stage reduced `fmax` from `0.619` to `0.309 eV/A`; after CI
  activation the residual eventually reached `0.213 eV/A` at step 30. The
  selected CI image had negative tangent curvature (`-4.954 eV/A^2`), proving
  the staged-CI branch executed, but the requested `0.05 eV/A` convergence was
  not reached.
- The path was monotonic and downhill (`Delta H=-0.099232 eV`), so the stored
  forward barrier is `0 eV` for this unconverged path and is not a physical
  transition-state result. This is a valid ABACUS/VCNEB execution diagnostic,
  not a completed material barrier.
- Records: `outputs/batio3_vcneb_pbe100_v2_summary.json`,
  `outputs/batio3_vcneb_pbe100_v2_summary.txt`, and
  `outputs/batio3_vcneb_pbe100_v2_manifest.json`.

## Core staged-CI and path-fold guard (`d9f3f82`)

- Added `run_vcneb(..., climb_after=N)`: the chain starts with ordinary NEB
  and enables CI after `N` completed optimizer steps.
- Extended `path_geometry_diagnostics()` with adjacent extended-coordinate
  segment lengths/cosines, zero-length segments, and optional fold rejection
  through `fold_cosine_threshold`.
- The full HF regression suite passes, including the staged-CI API and a
  folded-path negative test.  The default remains backward-compatible:
  diagnostics do not reject a path unless a fold threshold is requested.

## Random-path robustness and staged CI (`c0a6275`)

- The first direct-`climb=True` robustness run exposed a real failure mode:
  one deterministic perturbation converged to a folded, projection-stationary
  image chain with a spurious `0.7911 eV` barrier.  Its small NEB residual did
  not certify a physical MEP.
- The test was changed to the standard two-stage protocol: relax with
  `climb=False`, then restart on the same images with `climb=True`.  Four
  deterministic perturbation seeds and both `linear` and `log_strain` cell
  interpolation were run on HF using only `ToyPhaseTransition`.
- All 8/8 cases converged.  The largest barrier error was `1.2801e-4 eV`, the
  largest final generalized force was `1.9951e-3 eV/A`, and all final paths
  preserved the correct endpoint ordering.  The report also records raw
  atomic force, stress, cell-force, and final reaction-coordinate diagnostics.
- Record: `outputs/vcneb_robustness_staged_hf.json`.  The earlier direct-CI
  failure remains an important negative control and is not treated as a valid
  convergence result.

## Three mode-path semantics (`1badaa0`)

- Ran `examples/compare_mode_path_variants.py` on `hf` using only the coupled
  analytic `ToyPhaseTransition` calculator; no DFT job was submitted.
- With seven images and `fmax=0.002 eV/A`, the unconstrained, mode-guided then
  released, and strict coupled-mode branches give barriers `0.2500058`,
  `0.2500251`, and `0.2500000 eV`, respectively, against the exact `0.25 eV`
  reference.  All locate image 3 as the saddle with negative tangent
  curvature; the strict branch is already converged at initialization.
- Records: `outputs/mode_path_variants_hf.json` and
  `outputs/mode_path_variants_hf_manifest.json`.

## Algorithm-first periodic endpoint gauge validation (`d710cb9`)

- The first HF 7-image DFT trial was stopped after the initial force audit:
  the relaxed BaTiO3 endpoints had a periodic origin mismatch, and direct
  interpolation produced a shortest Ti-O distance of `1.328 A` with an
  initial generalized force of `39.018558 eV/A`.
- Implemented optional `align_translation=True`, which jointly optimizes a
  common periodic fractional translation and the element-grouped atom
  assignment.  The aligned calculator-free path has minimum distance
  `1.797 A` and maximum deformation `0.0631`; the CLI exposes hard geometry
  thresholds before any calculator is invoked.
- The complete regression suite passes, including the new periodic translation
  test.  The coupled analytic toy VCNEB/CI run recovers the known barrier as
  `0.250004 eV` versus `0.25 eV` and converges at step 47.
- The HF endpoint relaxations are retained as versioned fixtures and the
  subsequent DFT NEB run was intentionally cancelled until the algorithmic
  path checks are considered complete.  See
  `outputs/batio3/batio3_hf_algorithm_preflight.json`.

## Fixed-cell ASE reference comparison (`7259657`)

- The analytic coupled toy surface was run with the same seven images and
  climbing-image settings through ASE CINEB and VCNEB with `cell_mask=0`.
- ASE returned a barrier of `0.2811707727 eV`; VCNEB returned
  `0.2811614404 eV`, an absolute difference of `9.33e-6 eV`.  The reaction
  energy was `0.1774 eV` in both paths, and both identified image 4 as the
  saddle with negative tangent curvature.
- The machine-readable comparison is `outputs/fixed_cell_ase_comparison_hf.json`.

## Current-version analytic convergence matrix (`eb83831`)

- Re-ran the 54-case analytic matrix on `hf` using 5/7/9 images, springs
  `0.05/0.10/0.20 eV/A^2`, `cell_scale=4/5/6 A`, and FIRE/LBFGS, with
  `fmax=0.005 eV/A` and 300 optimizer steps.
- All 54 cases converged.  The known `0.25 eV` barrier was recovered with
  zero error at the stored output precision; the largest final generalized
  force was about `0.004988 eV/A`.
- Records: `outputs/vcneb_algorithm_matrix_hf.json` and
  `outputs/vcneb_algorithm_matrix_hf_manifest.json`.

The tangent regression now also exercises a non-collinear path in monotonic,
peak, and valley energy regimes.  It confirms the expected energy-weighted
improved-tangent direction in the same extended coordinate space used by the
spring test.

The strict-mode regression now includes a coupled atomic/cell mode whose
endpoint displacement lies exactly in the allowed subspace.  CI-VCNEB recovers
the analytic barrier and saddle residual without requiring an artificial
optimizer step when the initial path is already converged.

## HfO2 linear versus logarithmic-strain preflight (`403d06b`)

- Before the run, `cu17` was checked at `2026-09-12T08:15:51+08:00`:
  40 cores, load averages `0.08/0.03/0.08`, and no active DFT process.
- The calculator-free seven-image comparison used the mapped 12-atom HfO2
  fixture and held the endpoint order, MIC convention, and image count fixed.
  Both `linear` and `log_strain` paths were valid and retained positive cell
  determinants.
- The minimum interatomic distance changed only from `2.013034 A` (linear) to
  `2.012914 A` (logarithmic strain), while the maximum deformation norm stayed
  `0.049608`.  Thus changing cell interpolation alone does not resolve the
  current mechanism-path problem; no new DFT run was started on this result.
- The post-run check at `2026-09-12T08:17:28+08:00` found no VASP, ABACUS or
  MPI process remaining.

Compact records: `outputs/hfo2_initial_path_comparison_cluster.json` and
`outputs/hfo2_initial_path_comparison_manifest.json`.

## Endpoint atom mapping API (`266c8ea`)

- Added `infer_atom_mapping()` and `validate_atom_mapping()` for explicit or
  element-grouped automatic endpoint mapping.  The automatic path solves a
  finite square assignment problem per element using periodic Cartesian
  distances and returns an auditable per-atom displacement report.
- `interpolate_vcneb(mapping="auto")` applies the inferred final-atom
  permutation before interpolation; the default `None` remains identity for
  backward compatibility.  Large reconstructive transitions are documented as
  cases where an explicit chemical mapping is still required.
- Added regression coverage for swapped endpoint element order, auto mapping,
  explicit mismatch rejection, and preservation of the initial atom order.
- After checking `cu17` at `2026-09-12T08:10:17+08:00` (40 cores, load
  averages `0.08/0.03/0.11`, no active DFT process), the synchronized source
  passed the complete regression suite on the cluster.  The post-run check at
  `2026-09-12T08:11:47+08:00` found no VASP, ABACUS or MPI process remaining.

## Cell interpolation strategy selector (`97240fe`)

- Added an explicit `cell_interpolation` argument to `interpolate_vcneb()`.
  The default `linear` path remains backward compatible; `log_strain` uses a
  matrix-log/geometric path for symmetric-positive deformation gradients, and
  a user callback can supply a project-specific deformation path.
- The logarithmic strategy is intentionally conservative: it removes the
  relative rigid rotation when `align_cells=True`, rejects non-symmetric or
  non-positive deformation gradients, and leaves the existing positive-cell
  determinant guard in place for every image.
- Added regression coverage for exact endpoints, the geometric midpoint,
  callback execution, and rejection of an unaligned rigid rotation.  The full
  local regression suite passes.
- After checking `cu17` at `2026-09-12T08:01:47+08:00` (40 cores, load
  averages `0.08/0.03/0.16`, no active DFT process), the synchronized source
  passed compileall and the complete regression suite on the cluster.  The
  post-run check at `2026-09-12T08:03:00+08:00` found no VASP, ABACUS or MPI
  process remaining.

## Formal six-strain finite-difference validation (`404cbad`)

- Before the run, `cu17` was checked at `2026-09-12T07:48:54+08:00`:
  40 cores, load averages `0.00/0.02/0.32`, and no active DFT process.
- The latest source passed `compileall` and the complete regression suite on
  `cu17`, including the nonorthogonal cell-force regression.  The formal
  six-direction strain report then scanned three diagonal and three symmetric
  shear perturbations at six central-difference steps for zero pressure and
  `0.5 GPa`.
- The stable step was `epsilon=1e-6`.  The maximum six-component error was
  `3.1559e-11 eV` at zero pressure and `1.3529e-10 eV` at `0.5 GPa`.  The
  smallest step showed the expected round-off increase, so the report keeps
  the full step scan rather than presenting a single tuned number.
- The test uses a rotation-invariant metric cell energy.  This is important:
  ASE's symmetric Cauchy stress is appropriate for the physical symmetric
  strain directions, whereas a deliberately rotation-sensitive artificial
  energy would also contain antisymmetric work that the stress tensor cannot
  represent.
- The post-run check at `2026-09-12T07:51:03+08:00` found no VASP, ABACUS or
  MPI process remaining.

Compact records: `outputs/vcneb_finite_difference_report.json` and
`outputs/vcneb_finite_difference_manifest.json`.

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

## Current cluster regression (`799bc3e`)

- Before the run, `cu17` was checked at `2026-09-12T07:30:28+08:00`:
  40 cores, load averages `0.08/0.12/0.97`, and no active DFT process. The
  complete compile and non-DFT regression suite passed. A second check at
  `2026-09-12T07:31:23+08:00` found no VASP/ABACUS/MPI process remaining.
- The current suite includes the new path diagnostics, nonorthogonal mode
  projection, CI saddle, pressure, restart, parallel ownership, calculator
  contract, and ABACUS input-validator checks. The largest reported force
  transformation errors remain below `1e-9` in the test units.

Compact record: `outputs/vcneb_current_cluster_regression_manifest.json`.

## Latest cluster regression (`adaeb32`)

- Before the run, `cu17` was checked at `2026-09-12T07:37:08+08:00` and had
  40 cores with no active DFT process. The latest source passed compileall and
  the complete non-DFT suite, including path geometry preflight. A check at
  `2026-09-12T07:38:11+08:00` found no remaining VASP/ABACUS/MPI process.

Compact record: `outputs/vcneb_latest_cluster_regression_manifest.json`.

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
  to `https://github.com/xdzhu/varneb`.
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
