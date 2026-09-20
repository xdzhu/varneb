# Material-level multi-backend VARNEB production plan

This plan defines the two common benchmark paths requested for every backend:

| Case | Physical role | Canonical path | Production image count |
| --- | --- | --- | ---: |
| GaN B4 → B1 | finite-pressure path with an interior barrier | existing Qian endpoint mapping at 45.7 GPa | 29 total (27 interior) |
| BaTiO3 T → C | barrierless structural path | existing 100-Ry endpoint pair at 0 GPa | 7 total (5 interior) |

Both use ordinary VCNEB first, fixed endpoints, no CI in the initial
production run, and the project default convergence criterion `0.10 eV/Å`.
The GaN 29-image choice preserves the already validated Qian comparison
(27 interior images); the BTO 7-image choice is the established compact
barrierless path and leaves five independent workers.

## Common production contract

Each backend must pass, in order:

1. endpoint identity and atom-count gate;
2. fixed-endpoint static SCF/cutoff/k-point gate;
3. calculator-free path geometry preflight (`auto` mapping, MIC, translation
   alignment, `log_strain`, minimum-distance and deformation limits);
4. complete energy/force/stress contract for every interior image;
5. ordinary VCNEB convergence and per-image audit;
6. path plot, barrier/reaction enthalpy, and provenance record.

The controller evaluates each endpoint once and caches it. Only images
`1..n_images-2` are submitted to image workers. For the HF allocation, use
32 MPI ranks per image worker; `5 × 32 = 160` ranks for BTO and
`27 × 32 = 864` ranks for GaN are the upper bounds. If the allocation is
smaller, the executor runs the workers in waves without changing the path.

## Backend matrix

| Backend | ASE/VARNEB entry | GaN B4→B1 | BTO T→C | Promotion condition |
| --- | --- | --- | --- | --- |
| ABACUS | `vcneb.abacus` | run | existing production reference | ABACUS input and orbital/pseudopotential hashes |
| VASP | `vcneb.vasp` | existing accepted path; repeat under matrix manifest | run | VASP input contract and POTCAR hashes |
| QE | `vcneb.qe` / `make_ase_espresso_factory` | run with explicit Ga/N UPFs | run with explicit Ba/Ti/O UPFs | UPF approval plus cutoff/k-point convergence |
| CP2K | `make_ase_cp2k_factory` | run with GTH-PBE and matching basis | run with GTH-PBE and matching basis | basis/cutoff/SCF convergence |
| ABINIT | `make_ase_abinit_factory` | run with one reviewed PBE PSP family | run with one reviewed PBE PSP family | PSP format/XC and cutoff convergence |
| LAMMPS | `make_ase_lammps_factory` | only after a published GaN potential is selected | only after a published Ba-Ti-O potential is selected | potential/units/virial validation |

LAMMPS must not use the Ar Lennard-Jones smoke potential for either material.
If a suitable potential is unavailable on HF, that cell is explicitly blocked,
not replaced by a nonphysical surrogate. Energies from different calculators
must never be subtracted across rows; comparisons are shape/barrier trends
within a backend and against literature only after that backend's own gate.

## Execution and audit artifacts

The reusable driver is `examples/run_vcneb_ase.py`. Use a specialized factory
for QE/CP2K/ABINIT/LAMMPS and the existing VASP/ABACUS drivers where their
input contracts require it. Every production directory must retain:

```text
manifest.json
vcneb_preflight.json
initial-vcneb.traj
vcneb.traj
image_worker_manifest.jsonl
vcneb_summary.json
vcneb_barrier.png
images/image_0001/ ... image_<N-2>/
```

The material benchmark is promoted only when all rows that claim completion
contain complete per-image energies, forces, stresses, SCF status, endpoint
identity, input hashes, executable/module versions, and job IDs.
