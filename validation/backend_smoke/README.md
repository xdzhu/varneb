# Optional backend smoke evidence

These are contract-level checks on HF `hfacnormal01`, not material results.
Each job used one task, an isolated image directory, and the installed module
listed in the manifest.  A smoke passes only when ASE returns finite energy,
forces, and stress and the VARNEB capability preflight reports
`variable_cell=true`.

- LAMMPS: Lennard-Jones Ar2, `apps/lammps/23Jun2022-intel-21`, job
  `27740478`, success.
- CP2K: H2, GTH-PBE/DZVP-MOLOPT-SR-GTH, `apps/cp2k/intelmpi/2024.1-2021-plumed`,
  job `27740553`, success after the documented OT SCF setup.
- QE: H2 with the installed `H.SG15.PBE.UPF`, `apps/quantum-espresso/intelmpi/7.0`,
  job `27740644`, success with the static `scf` force/stress contract.

The H2 force is intentionally not interpreted physically; the case exists to
exercise the calculator contract.  For production use, a reviewed material
potential/basis, cutoff convergence, endpoint static gate, and a complete
VCNEB path are still required.
