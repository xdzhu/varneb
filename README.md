# VARNEB — variable-cell nudged elastic band

VARNEB is a calculator-agnostic Python toolkit for variable-cell NEB (VC-NEB)
paths between periodic crystal structures.  The public project name is
`varneb`; the import package remains `vcneb` for compatibility with the
original research scripts.

## Quick start

For an installed release, use `python -m pip install "varneb[plot]"`. From a
source checkout, this first example runs immediately and needs no DFT code or
material endpoint files:

```bash
python -m pip install -e ".[plot]"
varneb --version
python examples/run_toy_vcneb.py
```

It writes `toy_vcneb_run/` and should report a barrier near **0.25 eV**.
The potential is analytic: this verifies variable-cell forces and the run
interface, not a material prediction. To prepare a real material path next:

```bash
varneb backends
varneb optimizers
varneb init varneb.json
# edit initial, final, backend, workdir, and explicit calculator settings
varneb validate-config varneb.json
varneb prepare varneb.json
# inside a reviewed scheduler allocation, after loading the calculator:
varneb run varneb.json --execute
```

`prepare` only writes an initial chain and a geometry preflight report; it
does **not** run DFT. Set `calculator.kind` to `ase_class` or `factory`, pin its
`module:attribute` symbol, reviewed parameters and (for named DFT backends)
an explicit command before `run`. `run --execute` uses that same JSON and
starts calculators **in the current process**; it does not submit a scheduler
job or guess pseudopotentials/cutoffs. A seven-image path has five interior
workers; `image_workers: 0` runs them sequentially. Advanced restarts and mode
subspaces remain available via `python -m vcneb.material_runner` (the old
`examples/run_vcneb_ase.py` remains compatible). See
[`docs/USER_MANUAL.md`](https://github.com/xdzhu/varneb/blob/main/docs/USER_MANUAL.md) for the full contract. To run the
repository tests, install the additional `dev` extra and use `python -m pytest -q`.
`validate-config` checks JSON fields and types without reading endpoint files;
`prepare` is the endpoint and initial-geometry preflight. Use JSON booleans and
numbers, not quoted strings: typos and silent type coercions are rejected.

`varneb backends --json` prints the machine-readable backend capability table.
`varneb optimizers --json` prints the calculator-independent path strategy
table. The `backend` and `optimizer` fields in `varneb.json` are deliberately
orthogonal: changing FIRE/BlockFIRE/SplitFIRE/StagedFIRE/BFGS does not alter
the calculator profile, and changing the calculator does not silently alter
the path optimizer.
`varneb doctor` checks optional ASE adapters and executables visible in the
current shell; it does not launch a DFT calculation.  The generated
`varneb.json` is a deliberately small starting point, not a calculator
parameter guess.

## What is included

- `vcneb/` — the stable Python API, VC-NEB core, optimizers, provenance,
  modal/phonon analysis, and thin calculator adapters.
- `examples/` — calculator-free quickstarts and compact material fixtures.
- `docs/` — theory, backend contracts, validation plans, and the detailed
  user manual.
- `cluster/` — scheduler templates.  Production jobs run on HF through Slurm;
  endpoints are cached and only interior images are dispatched.
- `outputs/`, `validation/`, `benchmarks/` — auditable result and benchmark
  artifacts, never required for installing the library.

The old `README_VCNEB.md`, `run_NEB/`, and `run_VCNEB/` names remain as
compatibility entry points for existing research scripts.  New work should use
the public CLI, `vcneb` API, and the layout documented in
[`docs/REPOSITORY_LAYOUT.md`](https://github.com/xdzhu/varneb/blob/main/docs/REPOSITORY_LAYOUT.md).
Release artifacts and the deliberately tag-triggered PyPI workflow are
described in [`docs/RELEASE_PROCESS.md`](https://github.com/xdzhu/varneb/blob/main/docs/RELEASE_PROCESS.md).

## Minimal API

```python
from ase.io import read
from vcneb import interpolate_vcneb, run_vcneb

initial = read("initial/POSCAR")
final = read("final/POSCAR")
images = interpolate_vcneb(initial, final, n_images=7, mic=True)
# Attach a reviewed static calculator to each image first; for VCNEB it must
# provide energy, forces, and stress in ASE units and use a private directory.
chain, optimizer = run_vcneb(images, fmax=0.10, climb=False)
```

Start with ordinary NEB; enable a climbing-image refinement only after the
converged ordinary band has a genuine interior maximum above both endpoints.

`n_images=7` means **seven total images**: two fixed endpoints and five
interior images. The endpoints are not recalculated at every VCNEB update.
The default NEB force threshold is `0.10 eV/Å`. In a variable-cell run the
calculator must provide stress; with `cell_mask=0` the path is fixed-cell NEB
and needs only energy and forces. Cell updates remain owned by VARNEB, so QE
images must use `scf`, not `relax` or `vc-relax`. Endpoints are fixed and may
be reused from audited static calculations.

## Backends

ABACUS, VASP, QE, ABINIT, and CP2K have independently converged 45.7-GPa
GaN B4→B1 paths under their recorded first-principles contracts. ABACUS and
VASP additionally cover the largest set of material cases; CP2K also has a
converged barrierless BTO T→C case. LAMMPS is a separate classical-potential
adapter, not an interchangeable DFT validation. Every backend uses explicit
per-image directories and static force/stress preflight. See the
[`GaN CP2K case`](https://github.com/xdzhu/varneb/blob/main/examples/cases/gan_b4_b1_cp2k/README.md) for the complete
result and an important MPI-affinity warning before enabling multiple CP2K
image workers.

For HF module environments, inspect first and load only what the job needs:

```bash
ssh hf "module avail 2>&1 | grep -Ei 'lammps|quantum-espresso|cp2k|abinit|abacus|vasp'"
```

`varneb prepare` is calculator-free: it writes `initial-vcneb.traj` and a
`varneb_preflight.json` report before any DFT executable is called. See
[`docs/USER_MANUAL.md`](https://github.com/xdzhu/varneb/blob/main/docs/USER_MANUAL.md) for calculator-specific
factories, Slurm isolation, provenance, restarts, and modal analysis.

For any ASE calculator exposing energy and forces (plus stress for VCNEB), use
`examples/run_vcneb_ase.py`. Its `--cell-mode fixed` option enforces identical
image cells and does not request stress; the default `--cell-mode full` retains
variable-cell behavior. Specialized factories remain available when a code
needs a profile or legacy file protocol.

## Reproducible research figures

The [CPC manuscript package](paper/VARNEB_CPC/README.md) contains the
LaTeX sources, figure PDFs, plotted CSVs, and claim-to-evidence checklist.
Its two mode-surface examples have deliberately different scopes:

- BaTiO₃/ABACUS: a zero-pressure, `Q_y=0` symmetry-restricted variable-cell
  energy sheet with nine audited fitting nodes and two prospective DFT checks.
  It is not an unrestricted or finite-temperature potential surface.
- GaN/VASP: a 45.7-GPa, 600-eV **central** path-adapted enthalpy cut using
  18 path centers and 72 off-path static evaluations. The five-backend
  B4→B1 barrier paths separately validate the path controller; the
  two-coordinate cut does not certify a strict transition state.

Both contours display interpolation only inside their sampled domains;
barriers come from converged VCNEB chains, not from contour pixels.
