# BaTiO3 cubic to tetragonal fixture

This fixture contains a five-atom BaTiO3 endpoint pair for a compact variable-
cell phase-transition test:

- initial: cubic `Pm-3m`, source file `BTO-pm3m.vasp`;
- final: tetragonal `P4mm`, source file `BaTiO3.poscar`;
- endpoint composition/order: `Ba Ti O O O`;
- endpoint cells: cubic `3.977796 x 3.977796 x 3.977796 A` and tetragonal
  `3.990379 x 3.990379 x 4.102655 A`.

The source files were copied from the shared research filesystem:

```text
/home/zhuxd/deep/ABO3_data_gen/BTO/cubic_pm3m/BTO-pm3m.vasp
/home/zhuxd/deep/ABO3_data_gen/BTO/t_P4mm/BaTiO3.poscar
```

The same data set provides the ABACUS assets used by the smoke calculation:

```text
/home/zhuxd/deep/ABO3_data_gen/BTO/cubic_pm3m/STRU40.01x01x01/02.md
```

The direct identity-order path is rejected for production use because the
intermediate oxygen images nearly overlap (minimum distance about `0.072 A`).
The geometry-based mapping selects `[0, 1, 3, 2, 4]`, after which the seven-
image path has a minimum distance of about `1.347 A`.

This is a source fixture, not a claim that either endpoint is converged under
the final VCNEB production settings.  Endpoint relaxation, path geometry
preflight, and a reproducible ABACUS VCNEB run must be recorded separately.
