# VARNEB CPC submission gate

Current scope: one manuscript with no page limit and no separate supplement.
The 2026-09-28 build has eight main-text figures and one same-document
appendix figure. This gate supersedes the layout targets in the dated
`STAGE_GATE_2026-09-27.md`; it does not supersede raw-data audits.

## Claims the present evidence can support

| Manuscript claim | Evidence to check | Allowed conclusion and boundary |
| --- | --- | --- |
| Established NEB extended to variable cell through a backend-independent controller | `varneb_CPC.tex` Secs. 2–3; `tests/check_vcneb_forces.py`; fixed-cell and finite-difference regressions | A tested implementation and execution contract, **not** a new VCNEB physical formalism. |
| BTO/ABACUS barrierless T→C path and cubic Γ decomposition | Fig. 2/3 source data, `outputs/batio3_t_to_c_pbe100_dzp10au/` provenance, `bto_qy0_even_mode_two_holdout_gate_20260928.json` | Under PBE, 100 Ry, 10-au DZP, and a 4×4×4 electronic mesh, the seven-total-image path is monotonic; the 1×1×1 Γ basis explains atomic geometry. The 81-point frozen cut and nine-training/two-prospective-point `Q_y=0` restricted sheet are different objects. Neither is an unrestricted conditional PES, finite-temperature free energy, or activation barrier. |
| GaN 45.7-GPa B4→B1 multi-backend barrier topology | Fig. 4 source data, `evidence/gan_45p7_multibackend_vcneb_20260924.json`, `evidence/gan_45p7_full_image_input_contract_audit_20260929.md`, `benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json` | Five 29-total-image paths satisfy the common 0.10-eV/Å generalized-force criterion and share a dominant peak. The final visible inputs have consistent per-backend settings, but this is not a historical-call or raw-output audit. Cross-backend energies use separate endpoint baselines and calculator contracts. The 45.7-GPa pressure is adopted from Qian, not re-established as each PBE backend's coexistence pressure; Qian's PW91/ultrasoft/force-norm protocol differs. CP2K original-text atomic forces and final run-end markers are unavailable; cache-chain forces and raw energy/stress reconciliation do not replace them. |
| GaN atom–strain and mode interpretation | Figs. 5–8 source CSV/QA, `gan_600eV_joint_gamma_bridge_20260928.json`, `gan_600eV_ts_2d_9x9_refinement_20260928.json`, `gan_600eV_atomic_tube_dense_20260928.json` | At the original VASP/PBE 600-eV, 45.7-GPa contract, the local joint Hessian has a negative direction and sampled 81-point and central 90-point frozen cuts support a bounded mechanism interpretation. The energy–force/stress mismatch and incomplete stationary-point test preclude a strict variable-cell TS certificate. The central cut is not a whole-path or relaxed two-mode surface. |
| Case-specific convergence acceleration | Table 2, `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md`, Slurm launch audit | Matched-start ABACUS chains use the same first-crossing threshold of 0.10 eV/Å; BTO/HfO₂ counted-launch reductions are 78.5%/38.2%. Negative transfers and path differences remain visible. Do not advertise universal speedup or identical saddle basins. The pre-DFT feasibility gate is a separate safety mechanism; this ablation did not isolate its contribution to launch savings. |
| Installable and usable source package | `docs/RELEASE_NOTES_v0.0.2.md` local artifact smoke; `docs/USER_MANUAL.md`; `tests/test_cli_run.py`, `tests/test_config.py` | The clean local wheel/sdist build, fresh virtual-environment install, calculator-free seven-image `prepare`, and analytic toy run have passed. This does not certify the hosted PyPI artifact or a user's DFT executable. |

`MANUSCRIPT_EVIDENCE.md` and `FIGURE_LOGIC_AND_STYLE.md` are the current
figure-to-source maps. A pretty interpolated pixel is never an additional
DFT observation. The claim limits above must remain in the abstract, results,
captions, and conclusions when wording changes. The direct manuscript-number
regression `tests/test_cpc_numeric_claims.py` recomputes headline GaN barriers,
BTO blind errors and monotonic endpoint rise, Slurm launch reductions, and
the HfO$_2$ per-cell conversion from committed source tables.

## Remaining gates for this CPC submission

1. Replace author-order, affiliation, corresponding-author, CRediT,
   funding/computing-acknowledgement, and competing-interest placeholders
   with author-approved statements. Resolve the Program Summary's pending
   CPC Library field according to the submission workflow.
2. Obtain a final scientific read-through of every literature comparator,
   pressure and formula-unit normalization, figure caption, and Data
   availability statement. In particular, retain the BTO PBEsol distortion
   versus NEB distinction and the GaN/CP2K and TS limitations. The first
   DOI/source-page pass and five corrected author records are documented in
   `CITATION_AUDIT_2026-09-28.md`; this does not replace the final read-through.
3. On the final source revision, rerun the full test suite, manuscript claim
   check, LaTeX build, and rendered-page inspection. Record failures rather
   than treating an old green result as proof for changed files.

No new material calculation is a prerequisite for the **bounded claims
actually made in this draft**, unless a final source audit finds an error.
The unrestricted BTO lower envelope is obstructed at cubic C by the free
third soft Γ direction; strict GaN index-one TS certification, HfO₂ multi-mode
maps, bilayer hBN, and mode-guided release-and-refine remain follow-on
research, not retroactive claims of this CPC paper. Do not adjust BTO's
100-Ry/10-au-DZP contract or GaN's 600-eV contract to improve a plot.
