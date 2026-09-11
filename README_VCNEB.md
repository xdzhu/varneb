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

## Files

- `vcneb/core.py`: VC-NEB algorithm and optimizer-compatible object.
- `vcneb/modes.py`: mode-guided initial paths and modal path projections.
- `vcneb/vasp.py`: VASP input parsing and per-image calculator setup.
- `vcneb/abacus.py`: ABACUS calculator factory adapter.
- `examples/run_toy_vcneb.py`: analytic smoke test with a known 0.25 eV barrier.
- `examples/run_hfo2_t_po_model_vcneb.py`: mapped 12-atom HfO2 T -> PO geometry smoke test with a synthetic endpoint double-well calculator.
- `examples/run_vcneb_vasp.py`: VASP driver based on the existing endpoint layout.
- `examples/run_vcneb_abacus.py`: ABACUS driver skeleton.
- `scripts/setup_hfo2_t_po_validation.py`: builds the HfO2 T -> PO validation fixture from local source structures or portable copies.
- `scripts/validate_vcneb_inputs.py`: static dry-run validator for VASP/ABACUS VC-NEB image directories.
- `tests/check_vcneb_forces.py`: finite-difference checks for force/stress transforms.

## Mode-guided paths and component constraints

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

## Local checks

```bash
python tests/check_vcneb_forces.py
python examples/run_toy_vcneb.py
python scripts/setup_hfo2_t_po_validation.py
python examples/run_hfo2_t_po_model_vcneb.py
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

`examples/run_vcneb_abacus.py` expects an ASE ABACUS calculator.  On `cu05`,
both `ase.calculators.abacus` and `abacus_neb` are importable.  Before a real
run, fill material-specific settings in the `parameters` dict:

- `ecutwfc`
- `kpts`
- `pp`
- `basis`
- `pseudo_dir`
- `basis_dir`
- spin, smearing, van der Waals, and convergence settings

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

- Endpoints must have the same atom count, symbols, and atom order.
- The code removes a global cell rotation during interpolation, but production
  phase-transition work still needs careful endpoint matching.
- Stress can only drive physical strain components; arbitrary cell rotations
  are a gauge, not a real force degree of freedom.
- The current mode interface provides mode-guided seeds and post-processing;
  strict arbitrary mode-subspace MEPs still need an explicit constrained-MEP
  implementation and endpoint-subspace validation.
- This is a working prototype, not yet a published SSNEB implementation.  Treat
  DFT results as research data: compare against fixed-cell NEB and endpoint
  cell-relax results before trusting barriers.
