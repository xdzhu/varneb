# Variable-Cell NEB Prototype

This repository now contains a calculator-agnostic VC-NEB prototype in
`vcneb/`.  It is meant for crystal phase-transition barriers where the cell
changes along the path.

The project direction is a pure-Python `OpenVCNEB` toolkit: no MATLAB or USPEX
runtime is required.  VASP, ABACUS, and future calculators are external
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

The package is installable without MATLAB or USPEX:

```bash
python -m pip install .
vcneb --version
```

Before an optimization, `run_vcneb()` checks every image calculator for this
contract. Missing stress is a hard error because a variable-cell calculation
must not silently replace the cell force by zero. The same preflight report is
available through `inspect_calculator()` and `validate_image_calculators()`.

After a run, `chain.saddle_diagnostics()` reports the highest interior image,
its relative enthalpy, residual generalized force, tangent/normal force
components, and a local finite-difference tangent curvature.  The curvature is
only a path-local diagnostic, not a full Hessian or a proof of first-order
saddle character.

## Files

- `vcneb/core.py`: VC-NEB algorithm and optimizer-compatible object.
- `vcneb/modes.py`: mode-guided initial paths and modal path projections.
- `vcneb/vasp.py`: VASP input parsing and per-image calculator setup.
- `vcneb/abacus.py`: ABACUS calculator factory adapter.
- `vcneb/calculator.py`: capability preflight and image-aware calculator diagnostics.
- `examples/run_toy_vcneb.py`: analytic smoke test with a known 0.25 eV barrier.
- `examples/run_hfo2_t_po_model_vcneb.py`: mapped 12-atom HfO2 T -> PO geometry smoke test with a synthetic endpoint double-well calculator.
- `examples/compare_initial_cell_paths.py`: calculator-free comparison of linear and logarithmic-strain initial paths for any ASE-readable endpoint pair (defaults to HfO2); supports `--mapping auto`.
- `examples/run_vcneb_vasp.py`: VASP driver based on the existing endpoint layout.
- `examples/run_fixed_cell_ase_comparison.py`: ASE CINEB versus fixed-cell VCNEB comparison.
- `examples/run_vasp_single_image_smoke.py`: real VASP energy/force/stress smoke driver.
- `examples/run_vcneb_abacus.py`: ABACUS driver skeleton.
- `examples/relax_abacus_endpoint.py`: independent ABACUS endpoint relaxation
  with optional variable-cell filtering and explicit optimizer step control.
- `scripts/setup_hfo2_t_po_validation.py`: builds the HfO2 T -> PO validation fixture from local source structures or portable copies.
- `scripts/validate_vcneb_inputs.py`: static dry-run validator for VASP/ABACUS VC-NEB image directories.
- `tests/check_vcneb_forces.py`: finite-difference checks for force/stress transforms.

## Cluster execution policy

For long NEB or DFT runs, use the shared cluster only. Log in through `235`,
inspect the load and existing processes on a candidate node, and run only on a
currently idle node among `cu17`, `cu22`, `cu23`, `cu24`, `cu25`, and `cu26`.
Each node has 40 cores, and the nodes share the project directory and software
environment, so a source sync to the shared `/home/zhuxd` path is sufficient.
Do not run long jobs on `235` or in local WSL, and do not use `qsub`/PBS for
this project.

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

## Local checks

```bash
python tests/check_vcneb_forces.py
python examples/run_toy_vcneb.py
python scripts/setup_hfo2_t_po_validation.py
python examples/run_hfo2_t_po_model_vcneb.py
python examples/run_abacus_single_image_smoke.py
python examples/relax_abacus_endpoint.py --help
python examples/run_vasp_single_image_smoke.py --help
python examples/run_fixed_cell_ase_comparison.py
python examples/run_vcneb_convergence.py
python examples/run_release_and_refine.py
python examples/run_finite_difference_report.py
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
  --fmax 0.05 \
  --steps 300 \
  --vasp-bin /home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std \
  --ncores 40
```

The script reads `CONTCAR`/`POSCAR` endpoints, interpolates fractional
coordinates and cell deformation, copies `POTCAR`, and uses static single-point
VASP settings with `IBRION=-1`, `NSW=0`, `ISIF=2`, `ISYM=0`.

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

For example, the HfO2 Dojo-FR setup on the shared cluster uses:

```bash
python examples/run_vcneb_abacus.py \
  --initial validation/hfo2_t_to_po/image_00/POSCAR \
  --final validation/hfo2_t_to_po/image_06/POSCAR \
  --workdir validation/hfo2_t_to_po/abacus_vcneb_smoke \
  --command "mpirun -np 40 /home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus" \
  --pseudo-dir /home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential \
  --basis-dir /home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/selected_Orbs \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_7au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_7au_100Ry_2s2p1d.orb \
  --ecutwfc 60 --kpts 1 1 1
```

The driver enforces `cal_force=1`, `cal_stress=1`, and `out_stru=1`, which are
required for VC-NEB.

ABACUS runs support the same `--resume` and `--resume-trajectory` options as
the VASP driver.

For an ABACUS template directory, the dry static check is:

```bash
python scripts/validate_vcneb_inputs.py \
  --mode abacus \
  --template-dir path/to/abacus/template \
  --image-root path/to/vcneb/images
```

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
