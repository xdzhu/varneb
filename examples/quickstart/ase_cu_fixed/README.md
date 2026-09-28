# Complete external-code-free CLI smoke

Run from the repository root after installing VARNEB and ASE:

```bash
varneb validate-config examples/quickstart/ase_cu_fixed/varneb.json
varneb prepare examples/quickstart/ase_cu_fixed/varneb.json
varneb run examples/quickstart/ase_cu_fixed/varneb.json --execute
```

This three-total-image, fixed-cell Cu path uses ASE's inexpensive EMT
calculator. It exercises the same JSON, endpoint mapping, preparation,
calculator dispatch, fixed-endpoint cache, and summary path as a material run,
but it connects two translated representations of the same periodic state.
The expected zero barrier is **not** a physical phase-transition prediction.
No DFT program, pseudopotential, stress, cluster allocation, or scheduler
submission is needed. `--execute` remains required because the command has
the same safety boundary for all calculators.

The summary is `examples/quickstart/ase_cu_fixed/run/vcneb_summary.json`.
With the specified one-step budget it should report `status: converged` and
`final_max_generalized_force_eV_per_A` below the default 0.10-eV/Å gate.
A second `run` in the same workdir is refused to protect its records. For
another trial, change `workdir` in a copy of the JSON rather than overwriting
the first run.
