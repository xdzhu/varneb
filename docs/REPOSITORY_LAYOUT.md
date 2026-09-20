# VARNEB repository layout and migration policy

The repository is a research codebase with older production runs mixed into
the source tree.  A clean layout must not silently delete or relocate a
restart directory, so cleanup is staged and backwards compatible.

## Canonical public tree

```text
vcneb/                 Python package and public calculator adapters
examples/quickstart/   tiny calculator-free examples
examples/cases/        compact, documented material fixtures
docs/                  theory, user manual, backend and validation contracts
cluster/               HF/Slurm templates only
scripts/               bounded analysis and case-preparation utilities
tests/                 calculator-free regression and adapter-contract tests
benchmarks/            reproducible timing/convergence inputs
validation/            accepted gates and provenance manifests
outputs/               reviewed result artifacts and source data
paper/                 manuscript source and figure data
```

## Compatibility and generated areas

`README_VCNEB.md`, `run_NEB/`, and `run_VCNEB/` are retained for scripts from
the original project.  The top-level `initial_state/`, `final_state/`,
`work/`, `tmp/`, `build/`, `dist/`, and `_varneb_work/` directories are run or
build state, not library source.  They are ignored or documented rather than
globally deleted because a user may have an active restart there.

New calculations should use a run root containing:

```text
runs/<system>_<endpoint>_<backend>_<images>_<git-short-sha>/
  manifest.json
  endpoints/initial/  endpoints/final/
  images/image_0001/ ... image_<N-2>/
  trajectory/  reports/
```

Only `images/image_0001` through `image_<N-2>` are dispatched.  Endpoint
directories hold immutable static results and are reused by the manager.

## Naming and provenance rules

Use lowercase `snake_case` for Python modules and scripts, a public backend
name (`vasp`, `abacus`, `qe`, `lammps`, `cp2k`, `abinit`) in run paths, and a manifest
with the git commit, calculator command, parameter digest, input hashes,
image count, threshold, and scheduler job IDs.  Never infer a physical
potential, pseudopotential, or unit system from a filename.

## Cleanup sequence

1. Add or update a canonical entry under `examples/quickstart`, `examples/cases`,
   or `docs`; keep a compatibility wrapper when an old path is used by tests.
2. Verify calculator-free tests and input preflight.
3. Move only an explicitly inventoried, inactive result directory; record the
   old-to-new path in its manifest.
4. Remove obsolete wrappers only in a separate, reviewed change.

This policy prevents a “tidy” refactor from destroying an expensive restart.
