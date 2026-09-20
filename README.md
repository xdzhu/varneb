# VARNEB — variable-cell nudged elastic band

VARNEB is a calculator-agnostic Python toolkit for variable-cell NEB (VC-NEB)
paths between periodic crystal structures.  The public project name is
`varneb`; the import package remains `vcneb` for compatibility with the
original research scripts.

## Quick start (calculator-free)

```bash
python -m pip install -e .
varneb --version
varneb backends
varneb init varneb.json
python examples/run_toy_vcneb.py
python -m pytest -q
```

`varneb backends --json` prints the machine-readable backend capability table.
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
[`docs/REPOSITORY_LAYOUT.md`](docs/REPOSITORY_LAYOUT.md).

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

VASP and ABACUS have the repository's strongest production validation.  QE,
LAMMPS, CP2K, and ABINIT are supported through optional ASE adapters with explicit
per-image directories and static force/stress preflight.  The adapter status
is intentionally reported as `validated`, `adapter`, or `planned` rather than
claiming a material result that has not been run and audited.

For HF module environments, inspect first and load only what the job needs:

```bash
ssh hf "module avail 2>&1 | grep -Ei 'lammps|quantum-espresso|cp2k|abinit|abacus|vasp'"
```

See [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md) for calculator-specific
factories, Slurm isolation, provenance, restarts, and modal analysis.
