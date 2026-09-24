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
| Figure source data | `figures/vcneb_*_source_data.csv` | committed with manuscript; BTO T-to-C plotted directly per one formula unit (0.08713-eV endpoint rise), whereas the HfO2 four-formula-unit cell is divided by four |
| VASP material-path benchmark | `docs/GAN_QIAN_2013_REPLICATION_PLAN.md`, `docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md`, `outputs/neb_literature_benchmarks/manifest.json` | Hefei VASP/PBE paths are archived for GaN B4/B1 tetragonal and hexagonal, GaN B3/B1, CdSe cell mapping, and HfO2 controls; per-path CSV/figures record normalization and literature provenance |
| VASP VCA backend cases (BST50, PZT50) and PTO control | `docs/vasp_vca_validation.md`, `results/bst50_vcneb_final/`, `results/pto_vcneb_final/`, `results/pzt50_vca_vcneb_final/` | backend-validation examples; all are monotonic T-to-C paths, so they do not establish a VCA transition-state barrier or enter the current manuscript without a dedicated figure/claim decision |
| VASP GaN B3-to-B1 production path | `docs/GAN_QIAN_2013_REPLICATION_PLAN.md`, `scripts/analyze_gan_qian_path.py`, `outputs/neb_literature_benchmarks/gan_b3_to_b1_literature_comparison_source_data.csv` | converged at `fmax=0.09903 eV/A`; the 0.9553 eV/GaN barrier exceeds the 0.57 eV/GaN literature value and is retained as a mapping-mechanism negative control |
| VASP CdSe rock-salt-to-wurtzite paths | `docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md`, `examples/cdse_sheppard_2012/`, `outputs/neb_literature_benchmarks/cdse_rs_to_wz_cell_mapping_literature_comparison_source_data.csv` | cell-mapping branch converged at `fmax=0.09589 eV/A`; atomic-mapping branch remains outside the manuscript evidence set |
| GaN 45.7-GPa multi-backend path | `evidence/gan_45p7_multibackend_vcneb_20260924.json`, `evidence/gan_cp2k_45p7_chain_step_0052.traj`, `figures/gan_multibackend_validation_source_data.csv` | ABACUS/VASP/QE/ABINIT/CP2K converged at 0.3274/0.3385/0.3297/0.2924/0.2928 eV/GaN; all dominant peaks are image 15; CP2K's 29-image force/stress and geometry audit passes |
| CP2K BTO path | `validation/backend_smoke/cp2k_bto_k4_vcneb_20260923.json` | converged at `fmax=0.08183 eV/A`; 0.06080-eV/f.u. endpoint rise and no interior barrier; topology check only |
| Architecture and multi-backend figures | `scripts/plot_cpc_method_and_multibackend_figures.py`, `figures/varneb_architecture.*`, `figures/gan_multibackend_validation.*` | Python-only editable SVG/PDF plus 600-dpi PNG; all five GaN paths have committed source data |
| Figure argument and style contract | `FIGURE_LOGIC_AND_STYLE.md` | each panel has a claim, source, normalization, and risk boundary; aligned axes, readable labels, framed legends, and explicit pending data |
| Real-material strict-mode release-and-refine | no authoritative material run | explicitly future work |

The manuscript must not claim algorithmic novelty over Qian VC-NEB,
G-SSNEB/SSNEB, or FD-NEB. Its software contribution is an open, pure-Python,
calculator-agnostic, restartable and evidence-oriented implementation.
