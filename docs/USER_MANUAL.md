# VARNEB detailed user manual

This manual is the advanced companion to the short workflow in
[`README.md`](../README.md). It describes the stable public interfaces; the
case plans under this directory contain system-specific scientific settings.

## 1. Install and inspect

```bash
python -m pip install -e .
varneb --version
varneb backends --json > backend-capabilities.json
varneb doctor
```

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
```

`prepare` is the source-of-truth input gate: it resolves paths relative to the
JSON file, checks endpoint composition/order, creates the calculator-free
initial trajectory, and writes `varneb_preflight.json`. It uses the tested
`log_strain`/automatic-mapping/MIC/translation-alignment defaults and can
reject a configured minimum-distance or deformation threshold before any
external executable is launched.

The material production driver `examples/run_vcneb_ase.py` accepts either an
ASE class (`--calculator module:Class`) or a VARNEB factory
(`--factory module:function`). The latter is preferred for QE, CP2K, ABINIT,
and LAMMPS because pseudopotential/profile and file-protocol details remain
explicit in the factory arguments.

`n_images` includes both fixed endpoints. A seven-image path therefore has
five worker images. The default ordinary-NEB criterion is `0.10 eV/Å`; set a
different value explicitly when a study requires it. Endpoints are not
recalculated at every iteration and may be supplied from an audited static
cache. Use `interpolate_vcneb(..., mic=True, cell_interpolation="log_strain")`
for the tested variable-cell initialization and keep the mapping report in the
run manifest.

## 3. Calculator contract

Every image calculator must provide `get_potential_energy()`, `get_forces()`,
and `get_stress()`. Stress is mandatory because it becomes the cell part of
the generalized VCNEB force. Use:

```python
from vcneb import validate_image_calculators
validate_image_calculators(images, require_unique_directories=True)
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
profiles, and stress convention. VARNEB only consumes ASE's energy/forces/
stress interface and validates that every image has a private directory.

An image directory must be private to one worker. The controller can then
retry or resume one image without mixing calculator files from another image.

### VASP, ABACUS, and QE

Use the existing factories in `vcneb.vasp`, `vcneb.abacus`, and `vcneb.qe`.
VASP inputs are frozen by the input contract; ABACUS input generation keeps
the calculator-specific files in the image directory; QE images must use
`calculation='scf'`, `tstress=True`, and `tprnfor=True`. QE UPFs must be
explicitly PBE-marked and pinned by the approved manifest.

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
        "cutoff": 300,
        "xc": "PBE",
        "max_scf": 200,
        "inp": "&FORCE_EVAL\\n  &DFT\\n    &SCF\\n      &OT\\n      &END OT\\n    &END SCF\\n  &END DFT\\n&END FORCE_EVAL",
    },
    command="cp2k_shell.psmp",
)
```

CP2K rejects `PROJECT` paths longer than 80 characters. VARNEB therefore uses
a deterministic short per-image alias when a shared HF path is long and copies
`.inp`, `.out`, and `.pos` back to the image directory. A CP2K smoke pass is
still not a basis/cutoff convergence result; complete that gate for a material
before comparing barriers.

### ABINIT

ABINIT uses the ASE `AbinitProfile` and requires an explicit pseudopotential
directory. The factory does not download or infer pseudopotentials:

```python
from vcneb import make_ase_abinit_factory

factory = make_ase_abinit_factory(
    parameters={"ecut": 10, "toldfe": 1.0e-6, "pps": "psp8", "kpts": (1, 1, 1)},
    command="srun --exclusive --ntasks=1 abinit",
    pp_paths="/path/to/reviewed/abinit-psp8",
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

## 6. Auditing and promotion

The minimum acceptance record contains complete per-image energy, forces,
stress, SCF status, endpoint identity/atom count, input and code hashes, path
mechanism, units, and the threshold used for convergence. Only after these
checks pass should a path be compared with literature or promoted into the
paper's benchmark tables. A monotonic path is reported as barrierless; do not
invent a climbing-image saddle for it.
