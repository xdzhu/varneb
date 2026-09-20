# VARNEB — Variable-Cell NEB (compatibility reference)

> The canonical user-facing README is [`README.md`](README.md).  This file is
> retained because older research scripts and links refer to `README_VCNEB.md`.

> VARiable-cell Nudged Elastic Band code with universal first-principles calculators

This repository now contains the calculator-agnostic VARNEB (VC-NEB) toolkit in
`vcneb/`.  It is meant for crystal phase-transition barriers where the cell
changes along the path.

The project direction is a pure-Python `OpenVCNEB` toolkit: no MATLAB or USPEX
runtime is required.  VASP, ABACUS, Quantum ESPRESSO (QE), and future calculators are external
backends behind an ASE-compatible energy/force/stress contract.

## ASE status

ASE's built-in `ase.mep.NEB` does not optimize variable cells in periodic
directions.  Its interpolation helper can interpolate cells, and ASE's cell
filters show the stress-to-generalized-force pattern, but the NEB object itself
rejects periodic images with different cells.  This prototype fills that gap by
optimizing an extended coordinate vector:

- atomic block: fractional coordinates mapped through the reference cell;
- cell block: deformation gradient `F`, where `cell = cell0 @ F.T`;
- force block: Cartesian forces and stress transformed into those coordinates.

The implementation works with normal ASE optimizers (`FIRE`, `BFGS`, `LBFGS`)
and with any ASE calculator that provides `energy`, `forces`, and `stress`.
For ASE's explicit line-search optimizer, use `optimizer="BFGSLineSearch"`;
set `line_search_retries=N` to retry only a reported `LineSearch failed!` at
the last complete chain state, shrinking the step cap by
`line_search_retry_factor` on each bounded retry.  Calculator failures and
other optimizer errors are never silently retried.

The package is installable without MATLAB or USPEX:

```bash
python -m pip install .
varneb --version
# ``vcneb`` remains a compatible legacy console alias.
```

Before an optimization, `run_vcneb()` checks every image calculator for this
contract. Missing stress is a hard error because a variable-cell calculation
must not silently replace the cell force by zero. The same preflight report is
available through `inspect_calculator()` and `validate_image_calculators()`.

After a run, `chain.saddle_diagnostics()` reports the highest interior image,
its relative enthalpy, residual generalized force, tangent/normal force
components, a local finite-difference tangent curvature, and whether an
interior peak rises above both endpoints.  A negative curvature alone is not
enough: a monotonic band can make the highest interior image look like a CI
candidate even though no transition-state barrier is present.  The curvature
and peak flags are path-local diagnostics, not a full Hessian or proof of
first-order saddle character.

## Files

- `[DFT queued] examples/cdse_sheppard_2012/`: independent 8-atom CdSe
  rock-salt -> wurtzite mapping cases on hf. See
  `docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md` for provenance and the
  PBE/PW91, empirical/DFT and small-barrier accuracy distinctions;
  submission is not a converged validation result.
- `vcneb/core.py`: VC-NEB algorithm and optimizer-compatible object.
- `vcneb/modes.py`: mode-guided initial paths and modal path projections.
- `vcneb/vasp.py`: VASP input parsing, exact-POTCAR per-image setup, and a
  one-virtual-site VCA adapter that pulls coincident component forces back to
  the physical site. See `docs/vasp_vca_validation.md` for validated systems
  and its deliberately narrow applicability boundary.
  The adapter writes lattice components with 17 significant decimal digits and
  verifies exact binary64 round-trip before VASP; the optional native lattice
  preflight uses that same serialization policy. This prevents last-digit ASE
  POSCAR rounding from changing VASP's Bravais classification without changing
  the manager geometry, symmetry tolerance, or physical path.
- `vcneb/abacus.py`: ABACUS calculator factory adapter.
- `vcneb/qe.py`: QE `pw.x` static-image factory; it rejects QE `relax` and
  `vc-relax` so cell updates remain manager-owned.
- `vcneb/phonons.py`: calculator-free Γ-point force-constant diagonalization
  and aligned-path normal-coordinate projections; it does not label a
  non-stationary NEB image as a phonon calculation.
- `vcneb/calculator.py`: capability preflight and image-aware calculator diagnostics.
- `vcneb/executor.py`: optional image-level concurrent calculator executor; the
  controller remains single-process and each external calculator job step must
  use isolated directories (and `srun --exclusive` on Slurm).
- `[unit] examples/run_toy_vcneb.py`: analytic smoke test with a known 0.25 eV barrier.
- `[model] examples/run_hfo2_t_po_model_vcneb.py`: mapped 12-atom HfO2 T -> PO geometry smoke test with a synthetic endpoint double-well calculator.
- `[model] examples/compare_initial_cell_paths.py`: calculator-free comparison of linear and logarithmic-strain initial paths for any ASE-readable endpoint pair (defaults to HfO2); supports `--mapping auto`.
- `examples/mode_template.json`: copy-and-edit JSON template for an atomic mode plus an optional cell deformation mode.
- `[production-template] examples/run_vcneb_vasp.py`: VASP driver based on the existing endpoint layout.
- `examples/analyze_path_gamma_modes.py`: postprocess a completed chain against a
  stationary-endpoint Gamma force-constant archive; it reports atomic normal
  coordinates separately from the variable-cell degrees of freedom. Supply the
  recorded endpoint mapping with `--reference-permutation` when NEB reordered
  same-species atoms, and `--reference-translation` when the initial-path
  metadata records a non-integer endpoint gauge translation.
- `[production-template] examples/run_vcneb_qe.py`: QE `pw.x` driver with a
  no-DFT `--validate-only` preflight; it supports 7 total images and
  manager-controlled interior-image workers. The preflight verifies every
  selected UPF is an in-directory PBE file with matching element metadata and
  records its SHA256; this is a provenance gate, not a cutoff-convergence claim.
- `[model] examples/run_fixed_cell_ase_comparison.py`: ASE CINEB versus fixed-cell VCNEB comparison.
- `[DFT-smoke] examples/run_vasp_single_image_smoke.py`: real VASP energy/force/stress smoke driver.
- `[production-template] examples/run_vcneb_abacus.py`: ABACUS driver skeleton.
- `[production-template] examples/relax_abacus_native.py`: ABACUS-native atomic/cell endpoint
  relaxation (`calculation cell-relax`, `relax_method bfgs`); ASE is used only
  for structure conversion and final `CONTCAR` export.
- `scripts/audit_native_endpoint.py`: read-only endpoint gate for native
  summaries (return code, convergence flag, composition, force and stress).
- `scripts/promote_hfo2_endpoints.py`: validates both endpoint summaries and
  atomically publishes only passing `CONTCAR` files to the production
  `relaxed_T/` and `relaxed_PO/` directories.
- `[DFT-smoke] examples/relax_abacus_endpoint.py`: ASE-driven endpoint adapter.  It is a
  numerically equivalent fallback to native `cell-relax` when the latter has
  step-control trouble; pass `--stress-kbar` to require a force-and-stress gate.
- `scripts/setup_hfo2_t_po_validation.py`: builds the HfO2 T -> PO validation fixture from local source structures or portable copies.
- `scripts/validate_vcneb_inputs.py`: static dry-run validator for VASP/ABACUS VC-NEB image directories.
- `scripts/validate_vasp_vcneb_static.py`: verifies that reused VASP endpoint
  inputs become static `IBRION=-1`, `NSW=0`, `ISIF=2` image calls with an
  explicit, fixed symmetry policy.
- `scripts/audit_vcneb_result.py`: calculator-free audit of a completed summary,
  including generalized-force, barrier, volume, geometry and interior-barrier gates.
  Use `--max-stress-kbar` when the production endpoint/path policy requires a
  common stress threshold in addition to the generalized-force gate.  Use
  `--fmax-target` only for an explicit re-audit under a different reported
  threshold, for example a loose `0.10 eV/A` NEB criterion.
- `scripts/compare_vcneb_images.py`: calculator-free comparison of completed
  5/7/9-image summaries, including shared calculator settings and explicit
  handling of consistent barrierless paths.  Pass
  `--allow-duplicate-image-counts` when comparing same-resolution variants
  such as `linear` versus `log_strain`; image-count convergence keeps the
  default uniqueness gate.
- `scripts/export_vcneb_metrics.py`: export completed summary diagnostics to a
  per-image CSV containing reaction coordinate, enthalpy, cell lengths/angles,
  volume, stress and NEB force components for plotting or paper tables.
- `scripts/export_vcneb_structural_metrics.py`: export selected MIC key-pair
  distances and optional extended-space mode projections from the latest
  complete trajectory; `--derive-endpoint-mode` is explicitly a structural
  diagnostic, not a phonon mode.
- `scripts/plot_material_validation_figure.py`: create an editable SVG/PDF
  BTO/HfO2 material-validation figure plus long-form path and barrier source
  data from completed summaries only; it labels literature values as external
  references rather than calculation replicas.
- `docs/batio3_validation_protocol.md`: fixed BTO settings, image-count
  convergence gates, CI staging and recovery/archive requirements.
- `docs/material_validation_guide.md`: accepted BTO/HfO2 material conclusions,
  reporting boundaries, and the calculator-free command that regenerates the
  editable validation figure and source-data tables.
- `paper/VARNEB_CPC/`: the CPC manuscript source, bibliography, figure PDF and
  figure source data; its evidence checklist prevents calculator-smoke or
  energy-scale comparisons from being presented as production benchmarks.
- `docs/vasp_vca_validation.md`: VASP VCA input/force contract, committed
  BST50/PTO/PZT50 summaries, and the boundary between isovalent VCA and
  aliovalent defect chemistry.
- `tests/check_vcneb_forces.py`: finite-difference checks for force/stress transforms.

## Cluster execution policy

For long NEB or DFT runs, use the shared cluster only. The current production
entry point is Hefei through `ssh hf` and Slurm. Before every submission,
inspect `sinfo`, `squeue`, and the account's own jobs, then select a suitable
partition, node allocation, and task count from the actual calculator cost.
Do not mechanically request 40 cores: small cells and cheap settings may use
fewer tasks, while a large k-point/plane-wave calculation may justify more.
Record the partition, node, task count, elapsed time, and input version with
each result. Do not run long jobs on the gateway or in local WSL, and do not
use `qsub`/PBS for this project.

## Mode-guided paths and component constraints

Before attaching an expensive calculator, an interpolated path can be checked
for cell and atom geometry:

```python
from vcneb import interpolate_vcneb, path_geometry_diagnostics

images = interpolate_vcneb(
    initial, final, n_images=9, mic=True,
    minimum_distance=1.0,
    maximum_deformation=0.8,
)
report = path_geometry_diagnostics(images)
```

The optional thresholds on `interpolate_vcneb()` raise a `ValueError` with the
image index when a path violates them.  The report itself is calculator-free
and contains the per-image volume, shortest periodic interatomic distance, and
deformation norm.

The cell path can use the historical linear deformation interpolation, a
logarithmic strain path for aligned symmetric-positive deformation gradients,
or a project-specific callback:

```python
images = interpolate_vcneb(
    initial, final, n_images=9, mic=True,
    cell_interpolation="log_strain",
)

def interpolate_cell(lam, deform0, deform1):
    return (1.0 - lam) * deform0 + lam * deform1

images = interpolate_vcneb(
    initial, final, n_images=9,
    cell_interpolation=interpolate_cell,
)
```

`log_strain` removes the relative rigid rotation when `align_cells=True` and
rejects unsuitable non-positive or non-symmetric deformation gradients.  A
custom callback must return a finite 3x3 deformation matrix; every resulting
cell is still checked for a positive determinant.

The endpoint atom order can be kept explicitly or inferred by element and
periodic geometry.  The inferred permutation is returned by
`validate_atom_mapping()` and should be saved with the run manifest:

```python
from vcneb import validate_atom_mapping

mapping_report = validate_atom_mapping(initial, final, "auto", mic=True)
images = interpolate_vcneb(
    initial, final, n_images=9, mapping=mapping_report["mapping"], mic=True,
)
```

Automatic mapping is a geometry heuristic, not a chemical identity proof.  For
large reconstructive transitions, inspect `mapping_report` or provide an
explicit final-atom permutation.

For periodic crystals whose endpoint coordinate origins differ by a lattice
translation, enable the joint translation/mapping alignment before any
calculator is attached:

```python
images = interpolate_vcneb(
    initial, final, n_images=9, mapping="auto", mic=True,
    align_translation=True, cell_interpolation="log_strain",
    minimum_distance=1.0,
)
```

This changes only the periodic coordinate gauge of the final endpoint.  The
translation is chosen together with the element-grouped assignment, because
independent MIC choices can otherwise create artificial short bonds in the
interior path.  Always inspect the resulting mapping and geometry report for
reconstructive transitions.

The mode feature is intentionally split into an initial-path generator and a
diagnostic projection.  The former bends the interior images with an
endpoint-zero envelope, then a normal unconstrained NEB calculation can relax
to the full-space MEP:

```python
from vcneb import Mode, mode_guided_path, project_path_onto_modes

mode = Mode.from_file("soft_mode.dat", n_atoms=len(initial))
images = mode_guided_path(
    initial,
    final,
    n_images=7,
    mode=mode,
    amplitude=0.20,
    envelope="sin",
    mic=True,
)
modal_coordinates = project_path_onto_modes(images, initial, mode)
```

`mode_guided_path()` accepts the same `cell_interpolation`, `mapping`,
`align_translation`, `minimum_distance`, and `maximum_deformation` controls as
`interpolate_vcneb()`, so a mode-guided path cannot silently bypass endpoint
mapping or calculator-free geometry checks.

For an explicit preflight before creating an optimizer:

```python
from vcneb import validate_image_calculators

reports = validate_image_calculators(images)
```

Each VASP or ABACUS image should have its own calculator directory. A custom
launcher can additionally call
`validate_image_calculators(images, require_directory=True,
require_unique_directories=True)`.

Text mode files contain `n_atoms` rows of three numbers or a flattened `3N`
vector.  JSON and NPZ files can also carry an optional `cell` 3x3 deformation
mode.  A normalized atomic mode has an amplitude in Angstrom; a cell mode is
dimensionless in deformation-gradient space.

For ABINIT-like fixed components, pass an `(n_atoms, 3)` mask.  `1` keeps a
component active and `0` keeps it fixed:

```python
chain = VCNEB(images, atom_mask=mask)
```

This is a constrained calculation when used during optimization.  A
mode-guided initial path by itself is not constrained and should be preferred
when the final result is intended to be a full MEP.

For strict mode-subspace dynamics, first build a basis in the same extended
coordinate space as VCNEB:

```python
from vcneb import Mode, build_mode_basis

basis = build_mode_basis(Mode([[1.0, 0.0, 0.0]]), initial)
chain = VCNEB(images, mode_basis=basis, constraint_mode="subspace")
```

`constraint_mode="subspace"` requires the endpoint displacement to lie in the
supplied basis and projects the initial interior images and every optimizer
update into that affine subspace.  Use `constraint_mode="projected"` when the
non-mode part of an existing initial path should remain fixed while only the
optimization update is projected.  `build_direction_basis()` provides the same
interface for one allowed direction per atom or arbitrary atomic direction
  combinations.  Call `direction_basis_conflicts()` before combining a
direction basis with component masks; it reports partially clipped columns,
fully inactive columns, and rank loss.  These are constrained transition
paths; a constrained saddle is not automatically a first-order saddle in the
full configuration space.

The three semantics are compared reproducibly by
`examples/compare_mode_path_variants.py`: an unconstrained path, a
mode-guided initial path followed by unconstrained VCNEB, and a strict
atomic-plus-cell mode subspace.  On the coupled analytic toy surface, using
seven images and `fmax=0.002 eV/A`, all three recover the same `0.25 eV`
barrier within `2.6e-5 eV`; the strict case is exactly on the known MEP.  The
HF record is `outputs/mode_path_variants_hf.json`.  This comparison is an
algorithm test: a strict constrained saddle must still be interpreted as a
saddle in the constrained space, not automatically as a full-space saddle.

The production ABACUS and VASP drivers expose the same mode controls.  Supply a
JSON/NPZ/text mode with `--mode`; add `--mode-guided` when the mode should bend
only the initial path, or select `--constraint-mode subspace`/`projected` to
project the VC-NEB dynamics.  For example:

```bash
python examples/run_vcneb_abacus.py \
  --initial relaxed_T/CONTCAR --final relaxed_PO/CONTCAR \
  --mode hfo2_endpoint_mode.json --mode-guided \
  --constraint-mode projected --n-images 7
```

`subspace` validates that the fixed endpoints are connectable by the supplied
mode basis; `projected` keeps the current interior path as the affine reference
and projects only subsequent updates.  The four Hefei Slurm templates accept
the equivalent `MODE_FILE`, `MODE_GUIDED`, `CONSTRAINT_MODE`,
`MODE_AMPLITUDE`, and `MODE_ENVELOPE` environment variables.  These options
are opt-in and do not alter the default unconstrained FIRE workflow.  Set
`VALIDATE_ONLY=1` on a Slurm template to exercise path/mode/calculator
preflight without launching ABACUS; this is useful before requesting a costly
production allocation.

For a perturbed or noisy initial path, use a two-stage CI protocol: first call
`run_vcneb(..., climb=False)` to relax the path, then call it again on the same
images with `climb=True` for saddle refinement.  Starting CI immediately can
select a wrong image and produce a folded, projection-stationary path.  The
robustness example tests four deterministic perturbation seeds with both
`linear` and `log_strain` cell interpolation; the staged protocol recovers the
analytic `0.25 eV` barrier in all eight cases.  Always inspect
`path_diagnostics()` in addition to the optimizer residual: the latter is the
projected NEB force and is not, by itself, a proof that the physical path is a
valid MEP.  The record is `outputs/vcneb_robustness_staged_hf.json`.

The same protocol is available directly through `run_vcneb(...,
climb_after=N)`: the first `N` completed optimizer steps use ordinary NEB and
CI is then enabled automatically.  `path_geometry_diagnostics()` reports
adjacent extended-coordinate segment cosines; pass
`fold_cosine_threshold=0.0` to `validate_path_geometry()` when a folded path
must be rejected.  The default remains diagnostic-only so existing workflows
are not silently changed.

If an explicit ASE line search is useful for a calculator, select
`optimizer="BFGSLineSearch"` and opt into a small bounded recovery budget:

```python
run_vcneb(
    images,
    optimizer="BFGSLineSearch",
    line_search_retries=2,
    line_search_retry_factor=0.5,
)
```

Only ASE's exact `LineSearch failed!` condition is retried.  Each retry starts
from the last complete chain state with smaller step caps; SCF, MPI, timeout,
invalid-cell, and unrelated optimizer errors remain visible to the caller.

For a physically interpretable CI result, require
`saddle_diagnostics()["has_interior_barrier"]` and inspect
`path_diagnostics()["ci_warning"]`.  If the band is monotonic or its interior
peak is below an endpoint, report a barrierless/unresolved path rather than
calling the highest interior image a transition state.

## Local checks

```bash
python -m pytest -q  # discovers and runs every non-DFT regression check
python tests/check_vcneb_forces.py
python examples/run_toy_vcneb.py
python scripts/setup_hfo2_t_po_validation.py
python examples/run_hfo2_t_po_model_vcneb.py
python examples/run_abacus_single_image_smoke.py
python examples/relax_abacus_endpoint.py --help
python examples/relax_abacus_native.py --help
python examples/run_vasp_single_image_smoke.py --help
python examples/run_fixed_cell_ase_comparison.py
python examples/run_vcneb_convergence.py
python examples/run_vcneb_robustness.py
python examples/run_release_and_refine.py
python examples/run_finite_difference_report.py
python scripts/audit_vcneb_result.py path/to/completed_vcneb_workdir
```

Expected toy output:

```text
barrier_eV=0.250004
delta_eV=0.000000
```

The HfO2 fixture validates a real 12-atom variable-cell path:

- initial endpoint: tetragonal `P4_2/nmc (137)`;
- final endpoint: polar orthorhombic `Pca2_1 (29)`;
- model VC-NEB smoke barrier: `0.800023 eV`.

These checks were copied to
`235:/home/zhuxd/abacus/agent-runs/20260618-vcneb` and passed on `cu05`.

## VASP run

From this repository root:

```bash
python examples/run_vcneb_vasp.py \
  --initial initial_state/relax \
  --final final_state/relax \
  --workdir run_VCNEB/run \
  --n-images 7 \
  --fmax 0.10 \
  --steps 300 \
  --vasp-bin /home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std \
  --ncores 8
```

The script reads `CONTCAR`/`POSCAR` endpoints, interpolates fractional
coordinates and cell deformation, copies `POTCAR`, and uses static single-point
VASP settings with `IBRION=-1`, `NSW=0`, `ISIF=2`, `ISYM=-1`, `SYMPREC=1e-4`.
The latter is a tested GaN regression policy, not a guarantee for every lattice
or VASP version. Overrides apply to the entire run, never as per-image repairs.
Each actual input write checks the frozen parameters/source fingerprints,
geometry and POSCAR roundtrip, and records `vasp_input_contract.json`.
See [the VASP input contract](docs/vasp_input_contract.md) for evidence and limits.

To build all seven isolated static image directories and record the effective
VASP parameters without launching VASP, append `--validate-only`. The report
is `vcneb_preflight.json`; `--image-workers N` subsequently enables a
manager-controlled pool for interior images when the launcher reserves
exclusive resources for each worker. `--image-workers 0` remains serial, but now
uses the same durable exact-state image cache and manifest as parallel runs.
The cache namespace automatically includes effective parameters and input hashes;
changing settings requires a new cache directory. Calculator-free preflight does
not certify VASP's Bravais initialization or SCF convergence.

To continue a stopped run from the latest complete chain snapshot:

```bash
python examples/run_vcneb_vasp.py \
  --initial initial_state/relax \
  --final final_state/relax \
  --workdir run_VCNEB/run \
  --n-images 7 \
  --resume
```

The trajectory stores one frame per image at each optimizer step.  Resume reads
only complete `n_images` frame groups and ignores an interrupted partial tail.
When appending to an existing run, snapshot files under `snapshots/` continue
from the next available `step_####` / `chain_step_####.traj` index.

Here `n_images` is the total chain length, including both fixed endpoints; for
example, `--n-images 7` means five optimizable interior images.

Before launching VASP, run a dry static check:

```bash
python scripts/validate_vcneb_inputs.py \
  --mode vasp \
  --template-dir initial_state/relax \
  --image-root validation/hfo2_t_to_po
```

## ABACUS run

`examples/run_vcneb_abacus.py` expects an ASE ABACUS calculator. The material
settings can be supplied on the command line, so the same VC-NEB driver can be
used with a different calculator configuration:

- `ecutwfc`
- `kpts`
- `pp`
- `basis`
- `pseudo_dir`
- `basis_dir`
- spin, smearing, van der Waals, and convergence settings

For a production HfO2 setup on the shared cluster, use the same 100-Ry,
full 10-au DZP policy as the archived high-accuracy paths:

```bash
python examples/run_vcneb_abacus.py \
  --initial validation/hfo2_t_to_po/image_00/POSCAR \
  --final validation/hfo2_t_to_po/image_06/POSCAR \
  --workdir validation/hfo2_t_to_po/abacus_vcneb_smoke \
  --command "srun -n ${SLURM_NTASKS:-8} /home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus" \
  --pseudo-dir /home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential \
  --basis-dir /home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_10au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc 100 --kpts 2 2 2 --scf-thr 1e-8
```

The driver enforces `cal_force=1`, `cal_stress=1`, and `out_stru=1`, which are
required for VC-NEB.

ABACUS runs support the same `--resume` and `--resume-trajectory` options as
the VASP driver. Before a costly run, `--validate-only` creates the image
directories, checks energy/force/stress capability and unique per-image
directories, and writes `vcneb_preflight.json` without launching ABACUS.
Completed runs include the calculator report, Slurm metadata and Git revision
in an atomically updated `vcneb_summary.json`; resuming does not overwrite the
original `initial-vcneb.traj`.

For image-level concurrency, pass `--image-workers N`.  The controller remains
one Python process and evaluates up to `N` independent interior calculators
concurrently; fixed endpoints are evaluated once and reused for the remainder
of that run.  Thus a seven-image run normally uses five workers, not seven.
on Slurm the calculator command must use `srun --exclusive`, and the allocation
must provide `N` times the MPI width requested by one calculator.  The default
is `0` (serial images), so ordinary runs retain the reference execution path.
`--image-retries K` retries only a failed image evaluation (default `0`); the
calculator must be restart-safe in its per-image directory when this is enabled.
Parallel runs append `image_worker_manifest.jsonl` (or the path supplied by
`--image-manifest`) after each controller evaluation batch, including status,
elapsed time and per-image attempt counts.
For recoveries that should reuse an exactly identical interior-image state,
pass `--image-cache-dir CACHE` and, when the calculator settings are shared,
`--image-cache-namespace NAME`.  Cache keys include the image index, species,
Cartesian positions, cell and PBC flags byte-for-byte; any coordinate or cell
change is a cache miss.  Entries are written atomically as compressed NumPy
records, the namespace is locked in `cache_metadata.json`, and the manifest
records `cache_hits`/`cache_misses`.  Use a separate cache directory (or a
distinct namespace) for different pseudopotentials, orbitals, cutoffs, k-meshes
or SCF settings; the cache never overrides the fixed-endpoint policy.
If one image fails in a concurrent batch, completed sibling images are still
written to the cache before the batch error is returned, so a resumed job does
not discard successful DFT work.
The bundled ABACUS adapter also disables ASE's process-global `ase_sort.dat`
write when the input atoms are already grouped by species (as in the BTO and
HfO₂ fixtures).  This avoids a thread-level parser race; non-identity atom
orders should use the serial image backend or an adapter that provides a
directory-local sort file.
For reproducible recovery from a known complete snapshot, combine
`--resume --resume-step N` with `--resume-trajectory`; negative `N` counts from
the end and `-1` means the latest complete chain.
The bundled VASP and ABACUS drivers write `vcneb_failure.json` in their work
directory on such an interruption.  Library callers can pass `failure_report=PATH` to `run_vcneb()`; if an optimizer
or calculator exception interrupts the run, an atomically written JSON report
records the exception, completed optimizer steps, the trajectory/snapshot
locations to use for recovery, and any existing calculator log paths used for
bounded failure classification.  If bounded line-search recovery was enabled,
the report also records each retry's reduced step cap.  The original exception
is still propagated when the retry budget is exhausted.
The Hefei BTO template `cluster/hf_batio3_vcneb_parallel.slurm` demonstrates
four 32-MPI workers in a 128-task allocation.

The HfO₂ production template `cluster/hf_hfo2_vcneb_distributed.slurm` uses
five isolated 32-MPI workers for seven total images (the two endpoints are
fixed and cached), with the 100-Ry and Orb-DZP-10au Hf/O inputs.

For an ABACUS template directory, the dry static check is:

```bash
python scripts/validate_vcneb_inputs.py \
  --mode abacus \
  --template-dir path/to/abacus/template \
  --image-root path/to/vcneb/images
```

To compare completed image-count branches without rerunning DFT, use
`scripts/compare_vcneb_images.py`.  The default barrier spread gate is 0.02 eV;
the discrete highest-image coordinate gate allows half of the smallest local
reaction-coordinate segment, since the sampled peak can move by a fraction of
one image when the band is refined.  Same-image-count parameterization or
optimizer variants must opt in to duplicate counts explicitly.

For convergence reporting, keep the run-time threshold and any later audit
threshold distinct.  The API and production-template default is
`fmax=0.10 eV/A`, the ordinary-NEB acceptance threshold for new VARNEB runs.
The archived HfO2 7/9-image and CI records were intentionally run at stricter
historical targets (`0.05` and `0.03 eV/A`, respectively); those values remain
provenance, not the current default.  For example, job 27687189 has
`final_max_generalized_force=0.086679 eV/A` and passes a re-audit with
`python scripts/audit_vcneb_result.py <workdir> --fmax-target 0.10`.  That
branch is therefore accepted at the default ordinary-NEB threshold, while its
different mechanism still keeps it out of the strict historical 7/9-image
convergence table.

The separately archived real-material subspace-to-release diagnostic
(`outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_mode_subspace_release_job27693085/`)
starts from a strict subspace chain and then releases every degree of freedom.
Its minimum was `0.097528 eV/A` at optimizer step 275, followed by 25 observed
steps of rebound to `0.108376 eV/A`.  Thus its saved step-275 chain may be used
as a `0.10 eV/A` ordinary-NEB snapshot; it is not a stricter historical
`0.05 eV/A` result or a substitute for the ordinary/CI image-count comparison.

## Current limitations

- Endpoints must have the same atom count and composition.  Identity order is
  supported for audited structures; `mapping="auto"` provides an
  element-grouped periodic geometry heuristic, and explicit permutations are
  recommended for reconstructive transitions.
- The code removes a global cell rotation during interpolation, but production
  phase-transition work still needs careful endpoint matching.
- Stress can only drive physical strain components; arbitrary cell rotations
  are a gauge, not a real force degree of freedom.
- Strict mode-subspace and projected-update constraints are available through
  `mode_basis`; direction-basis conflict diagnostics are available through
  `direction_basis_conflicts()`.  `VCNEB.path_diagnostics()` reports per-image
  cell metrics, raw atomic force/stress, cell force, and the true/spring/NEB
  force decomposition needed to audit a variable-cell path.  Generalized
  sparse projectors and release-then-refine workflows are still under
  development.  `examples/run_release_and_refine.py` demonstrates the
  supported two-stage workflow on an analytic coupled potential.
- This is a working prototype, not yet a published SSNEB implementation.  Treat
  DFT results as research data: compare against fixed-cell NEB and endpoint
  cell-relax results before trusting barriers.
