# VARNEB — variable-cell nudged elastic band

VARNEB is a calculator-agnostic Python toolkit for variable-cell NEB (VC-NEB)
paths between periodic crystal structures.  The public project name is
`varneb`; the import package remains `vcneb` for compatibility with the
original research scripts.

## Quick start (calculator-free)

For an installed release, start with `python -m pip install varneb`. To run
the bundled examples and tests from a source checkout, use:

```bash
python -m pip install -e .
varneb --version
varneb backends
varneb optimizers
varneb init varneb.json
# after editing the endpoint paths:
varneb validate-config varneb.json
varneb prepare varneb.json
python examples/run_toy_vcneb.py
python -m pytest -q
```

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
from vcneb import interpolate_vcneb, run_vcneb

images = interpolate_vcneb(initial, final, n_images=7, mic=True)
# Attach one static calculator with energy, forces, and stress to each image.
result = run_vcneb(images, fmax=0.10, climb_after=None)
```

The default NEB force threshold is `0.10 eV/Å`.  In a variable-cell run the
calculator must provide stress; cell updates remain owned by VARNEB, so QE
images must use `scf`, not `relax` or `vc-relax`.  Endpoints are fixed and may
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

For any ASE calculator exposing energy, forces, and stress, use
`examples/run_vcneb_ase.py`; specialized factories remain available when a
code needs a profile or legacy file protocol.
