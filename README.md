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
interface, not a material prediction. For the same public CLI used by material
runs, try the complete external-code-free
[`ASE/EMT example`](examples/quickstart/ase_cu_fixed/README.md), which exercises
`validate-config → prepare → run --execute` with bundled inputs.

To prepare a real material path:

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
job or guess pseudopotentials/cutoffs. It prints one result line and saves
the full `vcneb_summary.json`; use `--full-summary` to echo that JSON.
A seven-image path has five interior workers; `image_workers: 0` runs them
sequentially. Advanced restarts and mode
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
- `examples/` — external-code-free quickstarts and compact material fixtures.
- `docs/` — theory, backend contracts, validation plans, and the detailed
  user manual.
- `cluster/` — scheduler templates.  Production jobs run on HF through Slurm;
  endpoints are cached and only interior images are dispatched.
- `outputs/`, `validation/`, `benchmarks/` — auditable result and benchmark
  artifacts, never required for installing the library.

For a fixed substrate plane, use the explicit Python
[epitaxial boundary API](docs/EPITAXIAL_BOUNDARY.md): the endpoint BFGS filter
and NEB share the same open cell subspace. This is different from the
in-plane-variable slab/vacuum boundary and is not selected automatically by
the CLI.

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

In `varneb backends`, **validated** describes archived material evidence, not
readiness of the current machine. In particular, the clean ASE 3.29.0 wheel
used in our distribution smoke lacks `ase.calculators.abacus`; VARNEB does
not bundle that optional adapter or the ABACUS executable. Before an ABACUS
run, use a compatible ASE/ABACUS adapter in the job's Python environment or
provide an explicit custom calculator factory, check
`varneb doctor --backend abacus`, and pass a one-image energy/force/stress
preflight. An executable on `PATH` alone is not sufficient.

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
  energy sheet with 27 audited conditional-DFT nodes on a 9×3 grid. The
  original nine-node model predicted its 15 new nodes within 1.131 meV/BTO,
  below the declared 2-meV/BTO gate. It is not an unrestricted or
  finite-temperature potential surface. A
  separate frozen-cubic two-soft-mode figure in the main text contains a
  complete 17×17 grid of 289 *static DFT* points; the previous 9×9 grid
  predicted its 208 new nodes within 0.488 meV/BTO. Its smooth contour is
  only display interpolation, and the projected VCNEB path leaves that plane.
- GaN/VASP: a 45.7-GPa, 600-eV **central** path-adapted enthalpy cut using
  18 path centers and 144 off-path static evaluations (162 measured
  coordinates total). The 72 new interleaved nodes passed a 0.0073-meV/GaN
  prospective interpolation check. Two of them lie about 0.2 meV/GaN below
  their frozen-cell centerline references; this does not establish a lower
  relaxed MEP. The five-backend
  B4→B1 barrier paths separately validate the path controller; the
  two-coordinate cut does not certify a strict transition state. A separate
  local atom–strain cut near the highest image contains 289 measured static
  DFT enthalpies on a 17×17 grid, under the same 600-eV protocol. The prior
  9×9 interpolant predicted its 208 new DFT nodes within 0.00698 meV/GaN.

All plotted contours display interpolation only inside their sampled domains;
barriers come from converged VCNEB chains, not from contour pixels.

### Focused hafnia research track

A separate [HfO₂ research plan](docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md) tests
whether mode–strain control can reduce switching barriers without lowering
competing decay barriers from the same polar well. This is an ongoing,
bounded study, not a completed JCTC paper or a new-theory claim.
[Exact-SCF chain observations](benchmarks/hfo2_channels/20261008/chain_observations/README.md)
provide public energy/force/stress replay, continuous-reference mode
projections and independent cell-strain records without invoking DFT.
Their [representation/force audit](docs/HFO2_CONTINUOUS_CHAIN_OBSERVATIONS_2026-10-08.md)
documents why a common three-pattern plane does not fully represent the
current PO→M bottleneck. Provisional peaks remain distinct from converged
paths, certified saddles and independently validated energy predictions.
