# VARNEB CPC claim-to-evidence checklist

| Manuscript element | Source of truth | Current status |
| --- | --- | --- |
| Generalized atomic/cell coordinates and force transform | `docs/theory.md`, `tests/check_vcneb_forces.py` | regression covered |
| Fixed-cell reduction and CI semantics | `tests/check_vcneb_forces.py`, `examples/run_fixed_cell_ase_comparison.py` | regression covered |
| Restart, caching, and manager/worker ownership | `vcneb/executor.py`, `tests/check_vcneb_forces.py` | regression covered |
| BTO T-to-C 5/7/9 image comparison | `outputs/batio3_t_to_c_pbe100_dzp10au/bto_validation_provenance.json` | accepted production evidence |
| BTO cross-calculator endpoint identity reference | `outputs/batio3_t_to_c_pbe100_dzp10au/bto_abacus_endpoint_identity.json` | accepted ABACUS/PBE T-to-C endpoint fingerprints; required reference for QE/VASP gates |
| BTO cubic $\Gamma$ force constants, direct Phonopy eigenvectors, and T-to-C mode projection | `outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_phonon_provenance.json`, `bto_cubic_phonopy_gamma_eigenpairs.npz`, `bto_tetragonal_to_cubic_n7_gamma_modes_phonopy_direct.json` | 32-MPI job 27704219; ABACUS--Phonopy unit convention, complete-basis residual, and mapping/gauge are recorded |
| BTO $\Gamma$-mode mechanism figure | `scripts/plot_bto_gamma_mode_figure.py`, `figures/bto_gamma_mode_path_source_data.csv` | editable SVG/PDF plus PNG preview; single deterministic path, no statistical error bars |
| HfO2 T-to-PO ordinary/CI/control comparison | `outputs/hfo2_t_to_po_pbe100_dzp10au/hfo2_validation_provenance.json` | accepted production evidence |
| Figure source data | `figures/vcneb_*_source_data.csv` | committed with manuscript |
| VASP BTO T-to-C 7-total-image production path | remote cu17 run at commit `cdde89c`; pending committed summary/source data under `outputs/` | completed at `fmax=0.0951479 eV/A` with five interior images and cached endpoints; not manuscript-ready until repository provenance is synchronized |
| VASP VCA backend cases (BST50, PZT50) and PTO control | `docs/vasp_vca_validation.md`, `results/bst50_vcneb_final/`, `results/pto_vcneb_final/`, `results/pzt50_vca_vcneb_final/` | backend-validation examples; all are monotonic T-to-C paths, so they do not establish a VCA transition-state barrier or enter the current manuscript without a dedicated figure/claim decision |
| VASP GaN B3-to-B1 production path | `docs/GAN_QIAN_2013_REPLICATION_PLAN.md`, `scripts/analyze_gan_qian_path.py`, job 27721015 | numerically converged at `fmax=0.09903 eV/A`, but the 0.9553 eV/GaN barrier exceeds the 0.57 eV/GaN literature value; backend/restart evidence only pending mapping-mechanism diagnosis |
| VASP CdSe rock-salt-to-wurtzite paths | `docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md`, `examples/cdse_sheppard_2012/` | running two independent initial mappings; provisional snapshots are not manuscript evidence until the 0.10 eV/A path criterion and final audit pass |
| Real-material strict-mode release-and-refine | no authoritative material run | explicitly future work |

The manuscript must not claim algorithmic novelty over Qian VC-NEB,
G-SSNEB/SSNEB, or FD-NEB. Its software contribution is an open, pure-Python,
calculator-agnostic, restartable and evidence-oriented implementation.
