# GaN 45.7 GPa: compact final evaluated chains

These five 29-frame ASE trajectories are the **last complete chain block**
from each hashed, finished remote `vcneb.traj` named in
`../gan_45p7_remote_source_index_20260927.md`. The source trajectories
contained 45, 21, 50, 52, and 53 complete chain blocks for ABACUS, VASP,
QE, ABINIT, and CP2K, respectively. The extraction and cross-check are
reproducible with `scripts/audit_gan_multibackend_trajectory_curves.py`;
hashes and numerical diagnostics are in
`../gan_45p7_final_chain_audit_20260927.json`.

For every saved image, ASE energy, forces, stress and cell volume are finite.
At 45.7 GPa, `H=E+PV` from each chain agrees with its corresponding plotted
source CSV to less than `1.0e-12 eV/cell`. All five peaks are image 15 and
their barriers agree with the manuscript's compact evidence record.

This is a **trajectory-to-figure audit**. It does not independently parse
every original backend SCF/output log, certify the NEB projected-force norm,
or prove a common DFT total-energy zero across programs. Those raw-output
and executable/pseudopotential provenance gates remain open.
