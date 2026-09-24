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
   ASE ABINIT/LAMMPS/CP2K factories.  Existing VASP/ABACUS/QE modules remain
   thin, import-compatible adapters. `vcneb.config` provides the
   calculator-free JSON run contract and safe initial-path preparation.
   `make_ase_calculator_factory` is the generic escape hatch: any ASE
   calculator class/factory that returns energy, forces, and stress can be
   attached without adding a conditional to `vcneb.core`.
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
| QE | `vcneb.qe` | yes, `scf` only | approved PseudoDojo PBE GaN endpoints; 45.7 GPa production running |
| LAMMPS | `vcneb.backends` | yes via ASE | HF contract smoke; potential-specific validation pending |
| CP2K | `vcneb.backends` | yes via ASE | 800-Ry GaN endpoints and lazy 16-rank shell lifecycle validated; production running |
| ABINIT | `vcneb.backends` | yes via ASE | approved PseudoDojo PBE/PSP8 GaN endpoints; production running |

“Adapter” means the Python contract is implemented and tested; it does not
mean a particular potential, pseudopotential, cutoff, or literature barrier is
validated.  A case is promoted to material evidence only after complete SCF,
force/stress, endpoint identity, and provenance gates pass.

## User workflow

```text
varneb init varneb.json
edit endpoint paths and explicit calculator parameters
varneb validate-config varneb.json
varneb prepare varneb.json
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

To add a calculator, first try `make_ase_calculator_factory` with a prebuilt ASE
class/factory and an explicit profile/parameter dictionary.  Implement a
specialized factory only when the code's file protocol or input contract needs
it (as with ABINIT, CP2K, VASP, or QE).  Do not add calculator conditionals to
`vcneb.core`.  The factory must make the launch command, directory, unit
system, stress support, and input provenance explicit; the registry status
remains `adapter` until a real HF smoke and audit are committed.
