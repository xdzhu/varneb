# HfO2 T to PO VC-NEB fixture

This fixture uses a 12-atom sqrt(2) x sqrt(2) x 1 tetragonal HfO2 commensurate supercell and a mapped 12-atom polar/orthorhombic conventional HfO2 endpoint from the local ABACUS Born-charge examples. The tetragonal parent itself has 6 atoms; the 12-atom supercell is intentional so both VC-NEB endpoints have the same atom count and one-to-one mapping.

Generated files:
- `T_HfO2_12.vasp`
- `PO_HfO2_12_mapped.vasp`
- `initial_vcneb_path.traj`
- `image_00` ... `image_06` POSCAR files
- `mapping_report.json`

## Production VC-NEB record

The production convention is ABACUS/PBE with `ecutwfc = 100 Ry`, full 10-au
DZP Hf/O orbitals, a `2x2x2` k mesh and `scf_thr = 1e-8`. Both endpoints are
independently converged 12-atom Hf4O8 cells and are evaluated once/cached;
only interior images run under the distributed controller. The canonical
result record is
`outputs/hfo2_t_to_po_pbe100_dzp10au/hfo2_validation_provenance.json`.

| Path | Total / interior images | Job | Barrier (eV per 12-atom cell) | Reported generalized force (eV/A) |
| --- | --- | ---: | ---: | ---: |
| ordinary log-strain | 7 / 5 | 27678218 | 0.1567509 | 0.0455156 |
| ordinary log-strain | 9 / 7 | 27678407 | 0.1562092 | 0.0497191 |
| CI refinement | 7 / 5 | 27678507 | 0.1291722 | 0.0295475 |
| linear-cell control | 7 / 5 | 27678924 | 0.1596772 | 0.0495087 |
| subspace-to-release diagnostic (loose snapshot) | 7 / 5 | 27693085, step 275 | 0.1380114 | 0.097528 |

The ordinary 7/9 barrier spread is `0.0005417 eV`; the CI result is
`32.29 meV/f.u.` and should be compared only with an explicitly qualified
literature protocol. The endpoint-displacement mode-guided diagnostic is kept
separately because it satisfies only the loose `0.10 eV/A` criterion and is
not part of the strict ordinary/CI table.

The last row is also separate from that strict table.  It begins from the
real-material strict-subspace chain and releases all atomic-plus-cell degrees
of freedom before ordinary VC-NEB.  Its 300-step continuation reached a
minimum `0.097528 eV/A` at step 275, then was deliberately observed for 25
more steps and rebounded to `0.108376 eV/A`.  The archived step-275 chain is a
conventional loose-`0.10 eV/A` acceptance snapshot, not a `0.05 eV/A`
strict-convergence result.  It has a `0.1380114 eV` barrier per 12-atom cell;
the full trajectory, preflight and audit are under
`outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_mode_subspace_release_job27693085/`.

Use the material-validation figure and its source data under
`outputs/vcneb_material_validation_figure/` for the plotted comparison. It
states the LDA/QE/USPEX, 40-image literature convention separately from this
ABACUS/PBE record.
