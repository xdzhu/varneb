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
- ABINIT: H2 with the module-provided `H.psp8` test potential,
  `apps/abinit/intelmpi/8.6.1`, job `27740907`, success after loading the
  Intel runtime and bridging ABINIT 8.x's five-line `.files` protocol. The
  exact input, output, and pseudopotential hash are recorded in
  `abinit_hf_20260921.json`.

The H2 force is intentionally not interpreted physically; the case exists to
exercise the calculator contract.  For production use, a reviewed material
potential/basis, cutoff convergence, endpoint static gate, and a complete
VCNEB path are still required.

The original material-level CP2K endpoint gate is recorded in
`cp2k_bto_endpoint_gate_20260922.json`, but the subsequent path audit showed
that both endpoints had collapsed to the cubic phase.  It is retained as a
negative example: force/stress convergence alone does not establish phase
identity.  The corrected 4x4x4-k-point T/C endpoints pass both the explicit
phase gate and the current project gate (`fmax < 0.10 eV/A`, residual stress
`< 1.0 kbar`); see `cp2k_bto_k4_endpoint_gate_20260923.json`.

The corresponding original GaN endpoint gate is recorded in
`cp2k_gan_endpoint_gate_20260922.json`.  Job `27756555` later exposed a
provenance gap: translation alignment changed the effective CP2K final
endpoint, so the source-file static gate no longer described the structure
actually evaluated by the path.  Its effective endpoint has
`0.183 eV/A` force and about `6.04 kbar` stress and is rejected.  The exact
hash mismatch and repair are recorded in
`cp2k_gan_effective_endpoint_mismatch_20260923.json`.  Production drivers now
require static-summary hashes to match the mapped/aligned endpoints.

The first accepted high-pressure material path in this repair series is
ABACUS GaN B4-to-B1 job `27760826`: 29 total images, `45.7 GPa`, final
`fmax=0.096844 eV/A`, and barrier `0.327366 eV/GaN`.  The calculator contract,
endpoint residuals, path audit and literature comparison are recorded in
`abacus_gan_45p7_vcneb_20260923.json`.

An existing BaTiO3 cubic-to-tetragonal endpoint was also used for a QE
candidate static calculation (`27740755`) and a safe `auto`-mapping,
`log_strain`, minimum-distance preflight (`27740765`).  The detailed record is
`bto_qe_candidate_20260920.json`; it is deliberately not an approved QE
production benchmark because the Dojo UPF collection has not passed the
project's user-approval and cutoff-convergence gates.
