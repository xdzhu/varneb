# VARNEB detailed user manual

This manual is the advanced companion to the short workflow in
[`README.md`](../README.md). It describes the stable public interfaces; the
case plans under this directory contain system-specific scientific settings.

## 1. Install and inspect

```bash
python -m pip install -e ".[plot,mode]"
varneb --version
varneb backends --json > backend-capabilities.json
varneb optimizers --json > optimizer-capabilities.json
varneb doctor
```

The base package needs only ASE and NumPy; `plot` supplies Matplotlib for
figure-producing examples, and `mode` supplies SciPy for mode analysis. Add
the `dev` extra only when running the test suite. The README's analytic toy
run is the quickest executable check before configuring any DFT backend.

`doctor` only checks Python adapter imports and executables on `PATH`. On HF,
load a module inside the Slurm script, then run `varneb doctor` with the same
Python environment. Do not install a second copy of LAMMPS, QE, CP2K, or
ABINIT when the cluster module provides it.

## 2. Configuration and path semantics

Generate a starter configuration and validate it before attaching a
calculator:

```bash
varneb init varneb.json
# edit initial, final, backend, workdir, and explicit calculator parameters
varneb validate-config varneb.json
varneb prepare varneb.json
# from within a reviewed Slurm allocation, after loading the calculator module:
varneb run varneb.json --execute
```

`prepare` is the source-of-truth input gate: it resolves paths relative to the
JSON file, checks endpoint composition/order, creates the calculator-free
initial trajectory, and writes `varneb_preflight.json`. It uses the tested
`log_strain`/automatic-mapping/MIC/translation-alignment defaults and can
reject a configured minimum-distance or deformation threshold before any
external executable is launched.

`varneb run` uses the *same* configuration; `--execute` is required because it
can launch first-principles programs. It never submits a Slurm job by itself.
The initial path and preflight are checked again, and a prior execution in the
same workdir is refused rather than overwritten. For a standard ASE
calculator, set the following fields in `varneb.json` (shown with a
calculator-free toy class; substitute a reviewed material calculator and
settings for production):

```json
{
  "schema_version": 1,
  "backend": "ase",
  "initial": "initial/POSCAR",
  "final": "final/POSCAR",
  "workdir": "runs/example",
  "n_images": 7,
  "cell_mode": "fixed",
  "optimizer": "FIRE",
  "image_workers": 0,
  "calculator": {
    "kind": "ase_class",
    "symbol": "ase.calculators.emt:EMT",
    "parameters": {}
  }
}
```

For a named DFT backend (`abacus`, `vasp`, `qe`, `cp2k`, `abinit`, or `lammps`),
declare `calculator.command` explicitly and keep all scientific settings in
`calculator.parameters`; `kind: "factory"` and `calculator.factory_kwargs`
select a specialized image factory when profiles, potential paths, or launcher
protocols require it. The `backend` label and `optimizer` remain independent.
`cell_mode: "full"` requests stress and allows cell evolution;
`cell_mode: "fixed"` requires identical cells and zero pressure and requests
energy/forces only. Set `image_workers` to the number of concurrent *interior*
calculations; this must fit the scheduler allocation, including each worker's
MPI ranks. The initial two endpoints are evaluated once, not at every update.
The run manifest records the JSON hash, declared backend, calculator symbol,
parameters, and input geometry. Neither `validate-config` nor `prepare`
certifies that a DFT profile is physically suitable; the operator must review
pseudopotentials, cutoff, k mesh, stress units and the executable beforehand.

The packaged advanced driver `python -m vcneb.material_runner` accepts either an
ASE class (`--calculator module:Class`) or a VARNEB factory
(`--factory module:function`). The latter is preferred for QE, CP2K, ABINIT,
and LAMMPS because pseudopotential/profile and file-protocol details remain
explicit in the factory arguments. `examples/run_vcneb_ase.py` remains a
compatibility wrapper for existing source-checkout jobs, and advanced
`--resume-snapshot`/`--subspace-artifact` options remain on that driver.

`n_images` includes both fixed endpoints. A seven-image path therefore has
five worker images. The default ordinary-NEB criterion is `0.10 eV/Å`; set a
different value explicitly when a study requires it. Endpoints are not
recalculated at every iteration and may be supplied from an audited static
cache. Use `interpolate_vcneb(..., mic=True, cell_interpolation="log_strain")`
for the tested variable-cell initialization and keep the mapping report in the
run manifest.

The Python `VCNEB` and `run_vcneb` entry points default to ordinary NEB
(`climb=False`). Set `climb=True` only after an ordinary path has a genuine
interior energy/enthalpy peak above both endpoints. `climb_after` delays CI
when paired with `climb=True`; it is **not** an automatic peak or
transition-state certification gate.

The material command-line drivers also default to ordinary NEB. Use
`--climb` after checking the interior peak, or `--climb-after N` for an
explicit delayed CI run. `--no-climb` remains accepted for older job scripts
and overrides `--climb-after`; removing it alone does not enable CI.

## 3. Calculator contract

Every image calculator must provide `get_potential_energy()` and `get_forces()`.
Variable-cell paths additionally require `get_stress()`, which supplies the
cell part of the generalized VCNEB force. In fixed-cell mode, the cell is
identical in every image and VARNEB does not request or report stress. Use:

```python
from vcneb import validate_image_calculators
validate_image_calculators(images, require_unique_directories=True)
# For a fixed-cell path, pass require_stress=False, require_variable_cell=False.
```

For an ASE calculator not yet listed in the registry, use the generic adapter:

```python
from vcneb import make_ase_calculator_factory, attach_image_calculators

factory = make_ase_calculator_factory(
    MyAseCalculator,
    parameters={"cutoff": 400, "profile": reviewed_profile},
    command="srun --exclusive --ntasks=32 reviewed-code",
)
attach_image_calculators(images, workdir="runs/my_backend", factory=factory)
```

The calculator remains responsible for its own units, pseudopotentials,
profiles, and (when applicable) stress convention. VARNEB consumes ASE's
energy/forces/stress interface as required by the selected mechanical mode and
validates that every image has a private directory.

The generic driver defaults to full VCNEB. Select `--cell-mode fixed` for an
ordinary fixed-cell NEB; this rejects nonzero `--pressure-gpa`, unequal image
cells, and mode-subspace artifacts containing cell directions. A fixed-cell
calculator only needs energy and forces, including in threaded image workers.
This flag does not implement the symmetric in-plane strain coordinates needed
for a two-dimensional variable-cell comparison.

An image directory must be private to one worker. The controller can then
retry or resume one image without mixing calculator files from another image.

### Backend/optimizer separation

The production architecture has two independent selections:

1. The backend factory creates an isolated calculator for one image and must
   return energy and atomic forces; variable-cell runs also require cell stress.
2. The VARNEB controller applies the NEB tangent/spring projection and runs a
   calculator-independent optimizer (`FIRE`, `BlockFIRE`, `SplitFIRE`,
   `ImageScaledFIRE`, `StagedFIRE`, `BFGS`, `LBFGS`, or `BFGSLineSearch`).

Thus a convergence experiment changes only `optimizer` and its strategy
parameters; switching VASP, ABACUS, QE, CP2K, ABINIT, or LAMMPS changes only
the calculator profile/launcher. Backend-specific SCF retries are not path
optimizer retries, and neither is allowed to silently change the other.

### VASP, ABACUS, and QE

Use the existing factories in `vcneb.vasp`, `vcneb.abacus`, and `vcneb.qe`.
For the uniform ASE entry point, `vcneb.vasp.make_ase_vasp_factory` exposes
the same VASP input contract as an image factory, while
`vcneb.backends.make_ase_calculator_factory` can wrap any other ASE
calculator that returns energy, forces, and (for variable-cell use) stress. VASP inputs are frozen by
the input contract; ABACUS input generation keeps
the calculator-specific files in the image directory; QE images must use
`calculation='scf'`, `tstress=True`, and `tprnfor=True`. QE UPFs must be
explicitly PBE-marked and pinned by the approved manifest. The shipped GaN
PseudoDojo examples use `${VARNEB_PSEUDO_DIR}` rather than a user-specific
absolute path; set it to the reviewed UPF/PSP8 directory before launching QE
or ABINIT.

### LAMMPS

LAMMPS does not guess a potential. Supply `pair_style`, `pair_coeff`, units,
species order, and masses explicitly:

```python
from vcneb import make_ase_lammps_factory, attach_image_calculators

factory = make_ase_lammps_factory(
    parameters={
        "units": "metal",
        "specorder": ["Ar"],
        "masses": ["1 39.948"],
        "pair_style": "lj/cut 10.0",
        "pair_coeff": ["* * 0.0103 3.4 10.0"],
    },
    command="srun --exclusive -n 1 lmp_mpi",
)
attach_image_calculators(images, workdir="runs/lammps", factory=factory)
```

The factory sets an image-local `tmp_dir`. Do not use a shared LAMMPS
temporary directory in concurrent workers. The first HF contract smoke is
recorded in `validation/backend_smoke/hf_20260920.json`.

### CP2K

Use `make_ase_cp2k_factory` with an explicit shell command and basis/potential:

```python
from vcneb import make_ase_cp2k_factory

factory = make_ase_cp2k_factory(
    parameters={
        "basis_set": "DZVP-MOLOPT-SR-GTH",
        "pseudo_potential": "GTH-PBE",
        "cutoff_ry": 800,
        "xc": "PBE",
        "max_scf": 200,
        "inp": "&FORCE_EVAL\\n  &DFT\\n    &SCF\\n      &OT\\n      &END OT\\n    &END SCF\\n  &END DFT\\n&END FORCE_EVAL",
    },
    command="mpirun -np 16 cp2k_shell.psmp",
)
```

CP2K rejects `PROJECT` paths longer than 80 characters. VARNEB therefore uses
a deterministic short per-image alias when a shared HF path is long and copies
`.inp`, `.out`, and `.pos` back to the image directory. A CP2K smoke pass is
still not a basis/cutoff convergence result; complete that gate for a material
before comparing barriers.

`cutoff_ry` is converted explicitly at the backend boundary; a bare ASE
`cutoff` is in eV and must not be documented as Ry. The CP2K calculator is
lazy and ephemeral: attaching a 29-image path starts no shell, and each active
worker starts one MPI shell only for its own evaluation, then closes it. Thus
`IMAGE_WORKERS`, not total image count, bounds live MPI worlds. The shown
16-rank launch is the measured HF GaN per-image contract; re-benchmark it
elsewhere. **Do not equate an 80-task allocation with five isolated 16-rank
worlds.** A separate GaN continuation demonstrated that five plain `mpirun`
commands can pin all 80 processes to the same 16 CPU IDs, reducing each rank
to roughly 20% CPU and making a single SCF iteration take about 238 seconds.
For a portable safe starting point use one CP2K worker. Before enabling
multiple workers, run an affinity canary on the actual allocation and use a
launcher that gives each world a non-overlapping CPU set. The case-specific
record is in `examples/cases/gan_b4_b1_cp2k/README.md`. This is an execution
constraint, not an optimizer effect, and a run with overlapping CPU sets must
not enter an acceleration benchmark.

### ABINIT

ABINIT uses the ASE `AbinitProfile` and requires an explicit pseudopotential
directory. The factory does not download or infer pseudopotentials:

```python
from vcneb import make_ase_abinit_factory

factory = make_ase_abinit_factory(
    parameters={"ecut": 1400, "toldfe": 1.0e-7, "pps": "psp8", "kpts": (4, 4, 3)},
    command="mpiexec.hydra -bootstrap slurm -n 8 abinit",
    pp_paths="${VARNEB_PSEUDO_DIR}",
    pseudopotential_manifest="examples/material_profiles/abinit_gan_pseudodojo_sr_manifest.json",
)
```

The HF adapter smoke uses the module-provided `H.psp8` test potential. That
is an interface check only; a material calculation must pin the ABINIT
pseudopotential files, exchange-correlation functional, cutoff, k mesh, and
code revision before it can enter a cross-backend benchmark.

## 4. Modes and phonons

At a stationary endpoint, generate or load Gamma force constants and use
`examples/analyze_path_gamma_modes.py`. It reports:

- mass-weighted normal coordinates for every path image;
- degenerate mode subspaces rather than arbitrary vectors inside a degeneracy;
- segment tangent overlaps;
- VCNEB reaction coordinate and segment lengths;
- per-image squared-amplitude contribution fractions and dominant modes.

The cell degrees of freedom are reported separately. A Gamma mode projection
does not prove that a non-stationary path image is a saddle or that the mode is
the unique reaction coordinate. Preserve the endpoint mapping and any
fractional gauge translation in the command line and manifest.

### Reproduce a path projection

The archived BTO example uses a five-atom cubic reference, a `1×1×1` **phonon
supercell**, and a separate `4×4×4` electronic k mesh. Its force-constant
archive has the explicitly recorded `eV/angstrom.au` unit; do not use the
script's `eV/angstrom^2` default for that archive. After making the audited
endpoint structure and complete seven-image trajectory available locally,
the corresponding command is:

```bash
python examples/analyze_path_gamma_modes.py \
  --reference C_endpoint.cif --trajectory converged_vcneb.traj \
  --force-constants outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_force_constants.npz \
  --force-constant-unit eV/angstrom.au \
  --phonopy-eigenpairs outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_phonopy_gamma_eigenpairs.npz \
  --n-images 7 --reference-permutation 0,1,4,3,2 \
  --reference-translation 0.5,0.5,0.5605508507551892 \
  --output mode_projection.json
```

`C_endpoint.cif` and `converged_vcneb.traj` above are placeholders, not files
shipped with the wheel. The exact source job, original trajectory location,
atom permutation, and gauge shift for this BTO analysis are recorded in
[`bto_cubic_gamma_phonon_provenance.json`](../outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_phonon_provenance.json).
For another mapping, recompute those two alignment arguments rather than
copying the BTO values. The output amplitudes are mass-weighted atomic
projections in the fixed reference-cell convention; cell strain is not a
Gamma phonon amplitude.

### What a two-coordinate surface means

- A **frozen cut** evaluates actual static calculator points after varying
  only the declared coordinates. All omitted atomic modes and the cell follow
  a stated fixed rule. Interpolation may be drawn only within the sampled
  domain and must have an independent error check. Projecting a variable-cell
  path onto this plane does not put the path *on* the frozen surface.
- A **conditional surface** relaxes the explicitly declared remaining degrees
  of freedom at each fixed coordinate. Audit orthogonal force and stress,
  competing starts/branches, curvature where stability is claimed, and
  preselected holdout energies **and structures**. A converged local branch
  is not automatically the global lower envelope.
- A **path-adapted cut** uses path position `s` and a declared transverse
  coordinate. It can display a central segment without claiming to cover the
  endpoints or to use two global endpoint phonon modes. At nonzero external
  pressure compare `H = E + PV`, with one pressure and one formula-unit
  normalization throughout; the plotted 2D values are not finite-temperature
  free energies.

For BTO, cubic C has three independent unstable Gamma directions. Fixing only
`(Q_z,Q_x)` while freely minimizing the omitted `Q_y` therefore cannot yield
a two-dimensional lower envelope that passes through C at its cubic reference
energy. The manuscript's current BTO panel instead fixes `Q_y=0`, releases
the remaining atom--strain coordinates, and fits an even-mode model to nine
audited nodes. Two predictions fixed before independent ABACUS calculations
have absolute errors of 0.150 and 0.977 meV/BTO, within the declared
2-meV/BTO gate. The contour is drawn only inside the measured-coordinate
hull; it is a **symmetry-restricted local sheet**, not the global conditional
lower envelope. The older 59-point BTO figure is a distinct **frozen cubic**
cut; the five-point branch pilot demonstrates lowering when `Q_y` is released
but does not certify a continuous PES. Source tables and exact claim limits
are in [`FIGURE_LOGIC_AND_STYLE.md`](../paper/VARNEB_CPC/FIGURE_LOGIC_AND_STYLE.md).
For GaN, the audited 45.7-GPa surface is a **central frozen
atomic-transverse enthalpy cut** around a variable-cell path, not an
endpoint-spanning or orthogonally relaxed two-mode surface. It uses 18 path
centers and 72 off-path VASP/600-eV statics, with independently checked
along-path and inner-transverse interpolation errors of 0.166 and
0.159 meV/GaN. See
[`MANUSCRIPT_EVIDENCE.md`](../paper/VARNEB_CPC/MANUSCRIPT_EVIDENCE.md) for
the bounded point counts and holdout errors. Do not derive a new activation
barrier from either interpolated contour; obtain it from the converged VCNEB
image enthalpies under the same calculator contract.

An interior highest image is only a transition-state candidate. Calling it a
strict variable-cell TS additionally requires a sufficiently stationary
atomic-plus-cell enthalpy gradient, one negative direction of the **joint**
Hessian with converged finite differences, consistent energy/force/stress
derivatives, and two-sided descent to the intended basins. A fixed-cell
imaginary Gamma mode alone does not establish this. The current GaN evidence
supports local atom–strain coupling but is deliberately not labeled a strict
TS certificate.

## 5. HF/Slurm execution

Before a real job:

```bash
ssh hf "sinfo -p hfacnormal01"
ssh hf "squeue -u iai806"
ssh hf "module avail 2>&1 | grep -Ei 'lammps|quantum-espresso|cp2k|abinit|abacus|vasp'"
```

Use the cluster templates in `cluster/`. A single manager may launch exclusive
image-level steps, but `IMAGE_WORKERS * IMAGE_MPI` must not exceed the
allocation. Record job IDs, actual task counts, modules, input hashes, and the
VARNEB git revision. Do not submit a material path from a login shell, and do
not label a job as converged merely because it was submitted or because one
image returned a force.

`--validate-only` checks path geometry, endpoint hashes, profiles, and approved
manifests without instantiating calculators that can start external programs.
Runtime capability and private-directory checks run again inside the scheduled
production job before the first DFT evaluation.

## 6. Auditing and promotion

The minimum acceptance record contains complete per-image energy, forces,
stress, SCF status, endpoint identity/atom count, input and code hashes, path
mechanism, units, and the threshold used for convergence. Only after these
checks pass should a path be compared with literature or promoted into the
paper's benchmark tables. A monotonic path is reported as barrierless; do not
invent a climbing-image saddle for it.
