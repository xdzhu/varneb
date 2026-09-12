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

## Cluster diagnostic record

The shared ABACUS assets were used on the permitted 40-core nodes with PBE,
100 Ry, `1x1x1` k points, and SCF threshold `1e-6`.  Independent variable-cell
FIRE endpoint trials gave:

| endpoint | energy (eV) | volume (A^3) | max force (eV/A) | max stress (eV/A^3) |
| --- | ---: | ---: | ---: | ---: |
| cubic, 18 steps | -3721.389389 | 90.8105 | 0.0000 | 0.001249 |
| tetragonal, 40 steps | -3734.971506 | 89.2083 | 0.0562 | 0.006312 |

The subsequent seven-image `log_strain`, auto-mapped, no-climb VCNEB trial used
FIRE with `maxstep=0.02` for 40 steps.  It reached a maximum generalized force
of `3.3801 eV/A`, with image 1 the highest interior image but still
`11.7779 eV` below the cubic endpoint.  Thus the reported forward barrier is
zero for this unconverged low-cost setup and must not be interpreted as a
physical BaTiO3 phase-transition barrier.  The full machine-readable record is
`outputs/batio3/batio3_relaxed_endpoints_fire40_summary.json`, with the endpoint
summaries, optimizer logs, and initial-path identity/auto-mapping comparisons
in the same directory.

The original unrelaxed-endpoint trial is retained separately: LBFGS oscillated
after recovery, while FIRE reached a residual plateau near `1.25 eV/A`; this is
useful for testing restart and optimizer diagnostics, not for a scientific
barrier.

## HF algorithm-first preflight

The relaxed endpoints were recomputed on Hefei Slurm with 16 tasks per job,
without exclusive-node allocation.  The cubic endpoint converged in 6 FIRE
steps at `-3737.466765 eV` and `65.2741 A^3`; the tetragonal endpoint
converged in 39 steps at `-3737.566081 eV` and `67.6792 A^3`, with maximum
force `0.0041 eV/A`.

Before restarting a DFT band, the calculator-free path check was run with
`log_strain`, `mapping="auto"`, `align_translation=True`, and 7 images.  The
minimum periodic distance improved from `1.328 A` to `1.797 A`, with maximum
deformation `0.0631`; the configured hard thresholds are `1.6 A` and `0.10`.
The machine-readable record is
`outputs/batio3/batio3_hf_algorithm_preflight.json`.  The corresponding DFT
NEB job was intentionally cancelled before optimization while the algorithm
verification is being prioritized.
