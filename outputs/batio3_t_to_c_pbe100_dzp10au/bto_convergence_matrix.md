# BaTiO3 T→C ordinary VC-NEB convergence matrix

All runs use ABACUS PBE, 100 Ry, Ba/Ti/O `Orb-DZP-10au`, 4×4×4 k points,
SCF threshold 1e-8, `log_strain` cell interpolation, the same relaxed
endpoints, `maxstep=0.003`, and four concurrent 32-MPI image workers on
`hfacnormal01`.  CI was disabled for every ordinary-path run.

| total images | Slurm jobs | final max generalized force (eV/A) | reaction enthalpy (eV) | interior barrier | audit |
|---:|---|---:|---:|---|---|
| 5 | 27676251 | 0.0191050 | 0.0871289 | no | ok |
| 7 | 27676310 | 0.0196710 | 0.0871289 | no | ok |
| 9 | 27676513 | 0.0199074 | 0.0871289 | no | ok |

`scripts/compare_vcneb_images.py` reports `barrierless-consistent` for the
5/7/9 set: barrier spread is 0 eV and the highest-image coordinate spread is
within the protocol tolerance.  Increasing image count is therefore not
needed for this direct monotonic path.

## Reverse direction

The C→T 5-image ordinary path (`27676299` + `27676485`) converged to
`0.0195577 eV/A`, reaction enthalpy `-0.0872206 eV`, zero forward barrier, and
no interior barrier.

## CI gate

The ordinary paths have no interior energy maximum, so CI refinement is
explicitly withheld.  The decision and evidence are in `bto_ci_gate.json`.

Machine-readable endpoint/path settings, job IDs, resource model, reverse path,
CI gate and precision sensitivity are consolidated in
`bto_validation_provenance.json`.  This current 100-Ry record supersedes the
older diagnostic-only manifest under `outputs/batio3/`.
