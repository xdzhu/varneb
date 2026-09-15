# BaTiO3 T-to-C production VC-NEB case

This five-atom fixture is the compact perovskite validation case for VARNEB.
The physical reporting direction is tetragonal `P4mm` (T) to cubic `Pm-3m`
(C); the directory name is retained for compatibility with the original source
fixture. Both endpoints use the same ordered composition, `Ba Ti O O O`, and
the production path uses a geometry-audited, element-preserving mapping.

## Production convention

- ABACUS through ASE, PBE, `ecutwfc = 100 Ry`.
- Ba/Ti/O Dojo-NC-FR pseudopotentials and full `Orb-DZP-10au` orbitals.
- `4x4x4` k mesh, `scf_thr = 1e-8`, `scf_nmax = 150`, `log_strain` initial
  cell interpolation, FIRE with `maxstep = 0.003`.
- The independently relaxed endpoint structures are cached and fixed during
  each path run. Only the interior images are dispatched to image workers.
- Each production path uses four concurrent 32-MPI interior workers on
  `hfacnormal01`; a 7-total-image run therefore has five interior images but
  at most four active workers at once.

The canonical provenance is
`outputs/batio3_t_to_c_pbe100_dzp10au/bto_validation_provenance.json`; the
image-count table and CI decision are in the neighboring
`bto_convergence_matrix.md` and `bto_ci_gate.json` files.

## Accepted results

| Direction | Total images | Job | Final maximum generalized force (eV/A) | Reaction enthalpy (eV) | Interior barrier |
| --- | ---: | ---: | ---: | ---: | --- |
| T to C | 5 | 27676251 | 0.0191050 | 0.0871289 | no |
| T to C | 7 | 27676310 | 0.0196710 | 0.0871289 | no |
| T to C | 9 | 27676513 | 0.0199074 | 0.0871289 | no |
| C to T | 5 | 27676299 + 27676485 | 0.0195577 | -0.0872206 | no |

The three forward paths are consistently monotonic: the higher-energy cubic
endpoint is the maximum. Consequently, a climbing-image calculation is
correctly withheld rather than being used to manufacture a transition state.
The 5/7/9 comparison reports `barrierless-consistent` with zero barrier spread.

The `6x6x6`, `scf_thr = 1e-9` endpoint/path sensitivity check preserves this
barrierless topology but changes the reaction energy by about 15 meV. Report
the k-point convention with any absolute endpoint-energy comparison.

## Source structures and preparation

- `source_structures/BaTiO3_cubic.vasp` and
  `source_structures/BaTiO3_tetragonal.vasp` are portable source structures.
- `endpoint_relax_cubic_CONTCAR.vasp` and
  `endpoint_relax_tetragonal_CONTCAR.vasp` are the endpoint files used by the
  archived production chain.
- `hf_endpoint_relax_*` retains small endpoint-relaxation evidence, not a
  substitute for the canonical 100-Ry result record.

Before a fresh calculation, follow `docs/batio3_validation_protocol.md` and
run the static geometry preflight. Do not treat the early low-cost diagnostic
folders as physical barrier evidence; they remain only as restart and failure
handling fixtures.
