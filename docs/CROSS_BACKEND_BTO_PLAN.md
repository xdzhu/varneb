# BaTiO3 cross-backend validation gate

The optional backend smoke jobs prove the shared VARNEB contract, but they do
not establish that different calculators describe BaTiO3 consistently.  This
plan keeps the material case small and reuses the existing BTO endpoint/path
fixtures.

## Preconditions

1. **QE:** approve one explicit UPF collection, acquire the exact three files
   without renaming, verify MD5/SHA256, and run fixed-endpoint cutoff tests.
   The current SSSP 1.3.0 mixed-family register is intentionally a candidate,
   not authorization (`examples/bto_qe_sssp_1_3_pbe_precision_candidate.json`).
2. **LAMMPS:** select a published Ba-Ti-O potential that supplies forces and
   virial stress for the intended phase-transition physics.  A Lennard-Jones
   smoke cannot be used for BTO and VARNEB will not guess a potential.
3. **CP2K:** choose a reviewed GTH basis/potential set and converge cutoff and
   SCF settings on the existing BTO endpoint.  The CP2K short-project-path
   handling is already implemented and tested on HF.

## Execution order

```text
approved inputs -> endpoint static convergence -> endpoint identity gate
                 -> 7-image ordinary VCNEB (endpoints cached)
                 -> path/force/stress audit -> barrier and mechanism comparison
                 -> optional Gamma-mode projection
```

Use the same endpoint atom order, image count, threshold (`0.10 eV/Å`),
mapping, cell interpolation, and no-climb policy for every backend.  Record
calculator-specific cutoffs, k-point sampling, potential hashes, job IDs, and
the exact VARNEB revision.  Compare path shape and endpoint energies only after
each backend passes its own convergence gate; do not subtract energies from
different calculators as if they shared a zero.

## Promotion rule

The BTO cross-backend result can enter `benchmarks/` or the paper only when all
three backends have complete per-image SCF/force/stress records and an accepted
endpoint identity.  If a backend is unavailable or its potential is not
approved, retain its adapter smoke and mark the material comparison pending.
