# Quickstart examples

The smallest runnable example is the analytic model at
[`../run_toy_vcneb.py`](../run_toy_vcneb.py).  It needs no external
calculator and is the first check after installation.  The scripts at the
parent level are retained as compatibility entry points while new examples
are grouped here by intent.

For a complete `varneb validate-config → prepare → run --execute` walkthrough
without DFT or a scheduler, use [`ase_cu_fixed/`](ase_cu_fixed/README.md).

```bash
python examples/run_toy_vcneb.py
python examples/run_release_and_refine.py
python examples/compare_initial_cell_paths.py --help
```

First-principles material examples are separate from these local quickstarts;
they require reviewed inputs and the relevant external executables.
