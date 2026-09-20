# VARNEB architecture and user workflow

## Design boundary

The controller owns the path: endpoint mapping, cell interpolation, tangent
and spring decomposition, optimizer state, image cache, convergence, and
provenance.  A calculator owns one static evaluation of one structure and
returns only energy, Cartesian forces, and stress.  The controller never calls
an external executable directly.

```text
CLI/config -> case preparation -> endpoint cache
                         -> backend factory -> image workers
                         -> VCNEB controller -> trajectory/reports
                         -> modal/phonon analysis and plots
```

Every worker has an isolated directory.  This is essential for LAMMPS temporary
files, CP2K output prefixes, QE scratch files, and VASP POSCAR/INCAR/POTCAR
contracts.  The fixed endpoints are not dispatched during every iteration.

## Public layers

1. `vcneb.core`: calculator-independent VC-NEB state and optimization.
2. `vcneb.calculator`: capability checks and stable failure categories.
3. `vcneb.backends`: registry, backend status, common image attachment, and
   ASE LAMMPS/CP2K factories.  Existing VASP/ABACUS/QE modules remain thin,
   import-compatible adapters.
4. `vcneb.modes` and `vcneb.phonons`: mode-guided initial paths, strict mode
   subspaces, Gamma eigenvectors, path projection, degeneracy grouping, and
   tangent overlaps.
5. `examples/` and `scripts/`: user-facing drivers and bounded postprocessing;
   no calculator logic is copied into the core.

## Backend maturity contract

| Backend | Adapter | Variable-cell stress preflight | Material evidence |
|---|---|---:|---|
| VASP | `vcneb.vasp` | yes | production GaN/BTO/HfO2/CdSe evidence |
| ABACUS | `vcneb.abacus` | yes | production HfO2/BTO evidence |
| QE | `vcneb.qe` | yes, `scf` only | HF contract smoke; UPF/cutoff convergence pending |
| LAMMPS | `vcneb.backends` | yes via ASE | HF contract smoke; potential-specific validation pending |
| CP2K | `vcneb.backends` | yes via ASE | HF contract smoke; basis/cutoff validation pending |

“Adapter” means the Python contract is implemented and tested; it does not
mean a particular potential, pseudopotential, cutoff, or literature barrier is
validated.  A case is promoted to material evidence only after complete SCF,
force/stress, endpoint identity, and provenance gates pass.

## User workflow

```text
varneb init varneb.json
edit endpoint paths and explicit calculator parameters
varneb doctor --backend <name>
calculator-free geometry/preflight
static endpoint calculations and immutable cache
VCNEB with total image count (endpoints included)
audit convergence and interior barrier
project converged path onto endpoint Gamma modes
write figures, manifest, and reproducible report
```

The default force criterion is `0.10 eV/Å`; users can request a stricter
threshold, but a result is never silently re-labelled as converged under a
different threshold.  Climbing image remains an explicit second stage after
ordinary NEB and is skipped for a monotonic path without an interior barrier.

## Extension rule

To add a calculator, implement one factory and its preflight tests.  Do not add
calculator conditionals to `vcneb.core`.  The factory must make the launch
command, directory, unit system, stress support, and input provenance explicit;
the registry status remains `adapter` until a real HF smoke and audit are
committed.
