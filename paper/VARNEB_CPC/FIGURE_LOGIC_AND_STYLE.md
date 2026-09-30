# VARNEB CPC figure argument and evidence contract

This file describes the single two-column manuscript as of 2026-09-30:
**eight main-text figures, one appendix figure, and two main tables**. Figure
numbers below follow `varneb_CPC.tex`. Measured DFT markers, fitted contours, projected
paths, and external literature curves must remain visually and verbally distinct.

| Item | Claim and source | Boundary that the caption and text must retain |
| --- | --- | --- |
| Fig. 1, architecture | One manager owns the generalized path; workers, calculators, update policies, cache, and audits have separate interfaces. Source: `figures/varneb_architecture.pdf`. | An implementation diagram, not a measured scaling benchmark. |
| Fig. 2, BTO restricted sheet | ABACUS/PBE/100-Ry/10-au-DZP BTO at zero pressure. A measured 9×3 `Q_y=0` variable-cell sheet reuses 12 audited coordinates and adds 15 raw-audited conditional relaxations; the seven-image T→C path is projected onto it. Sources: `figures/bto_qy0_restricted_sheet_27_20260929_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/bto_qy0_27_node_sheet_20260929.json`. | Smooth color is shape-preserving interpolation of 27 measured nodes. The original nine-node even-mode model predicts the 15 new points within 1.131 meV/BTO, below its ±2-meV gate; two stress-gate misses required same-input optimizer continuation. Stable atomic modes and strain are released but not the plotted axes; the omitted cubic soft mode is fixed at `Q_y=0`. Neither an unrestricted conditional PES, finite-temperature FES, nor T→C activation barrier. |
| Fig. 3, BTO frozen two-soft-mode grid | The audited 9×9 grid is nested within a complete 17×17 grid using 208 new raw-audited ABACUS/PBE/100-Ry/10-au-DZP static calculations. Sources: `figures/bto_frozen_soft_mode_289_20260929_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/bto_transverse_soft_frozen289_20260929.json`. | The prior 9×9 interpolant predicts the 208 new nodes within 0.488 meV/BTO; the refined-grid interior LOO maximum is 1.61 meV/BTO and display minimum undershoot 0.637 meV/BTO. The projected VCNEB path leaves the frozen plane; this is distinct from Fig. 2's restricted sheet, not a relaxed PES or a barrier. |
| Fig. 4, GaN backends | Five converged 29-total-image B4→B1 paths at 45.7 GPa, with forward/reverse barriers referenced to their respective endpoints. Source: `figures/gan_multibackend_validation.pdf` and its committed source tables. The CP2K all-interior raw-energy/stress versus exact-cache report is `benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json`. | Energies are enthalpies per GaN. Literature curve is approximate digitization, not matched raw data; only the forward panel uses the published 0.34-eV/GaN guide. The pressure is adopted from Qian, whose QE/PW91/ultrasoft and RMS-force protocol is not the PBE/max-force contract here. Calculator contracts differ, so cross-code agreement is topological and approximate, not identity. CP2K raw atomic forces and final run-end markers remain unavailable, despite 27/27 energy/stress reconciliation. |
| Fig. 5, GaN local joint cut | Original VASP/PBE 600-eV joint atom–strain coordinates show a local negative and a positive curvature near the highest image. The final display interpolates 289 measured enthalpies on a complete 17×17 grid. Sources: `figures/gan_600eV_local_joint_dft289_20260930_v2_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_ts_2d_17x17_refinement_20260930.json`. | The prior 9×9 interpolant predicts the 208 independently measured new nodes with 0.00698-meV/GaN maximum and 0.00228-meV/GaN RMS error, passing predeclared 0.02/0.01-meV/GaN gates. The 25→81 check remains historical. This frozen local cut is neither an orthogonally relaxed or whole-path surface nor a strict TS certificate; the 0.0232-eV/Å energy–force mismatch remains. |
| Fig. 6, GaN Hessian and Γ projections | Frozen atomic, frozen strain, joint and strain-relaxed curvatures accompany the local Fig. 5 cut; endpoint optical-group projections describe the joint negative direction. Sources: `figures/gan_joint_mode_600eV_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json`. | The coordinate metric controls numerical curvature magnitudes. Endpoint optical overlaps are geometric, not mode-resolved barrier energies; the joint candidate is not a certified TS. |
| Fig. 7, whole-path GaN modes | Two independent endpoint Γ optical bases, B4-referenced strain, and reconstruction residuals across the mapped 29-image path. Sources: `figures/gan_gamma_path_600eV_source_data.csv` and `_qa.json`. | Five vectors in three leading optical groups capture >99.99% of this path's squared optical displacement; the B4/B1 bases are not mode-to-mode identified. Maximum residuals are 0.00043/0.00139 `sqrt(amu) Å`, while normal strains reach 34.0%, 24.3%, and 17.8%. This is geometry, not a decomposition of enthalpy. |
| Fig. 8, GaN central 2D cut | The 45.7-GPa central `H(s,q_perp)-H(s,0)` contour resolves the small transverse scale separately from the full path barrier in panel b. It uses 18 path centers and 144 atomic-only off-path statics under the byte-identical original VASP/PBE 600-eV contract. Sources: `figures/gan_600eV_atomic_transverse_162_20260929_v2_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json`. | The displayed 18×9 contour uses shape-preserving interpolation only inside the measured rectangle; the 72 interleaved new points pass a predeclared quartic-prediction gate with 0.0073-meV/GaN maximum error. Two measured nodes at image 19/20, +0.0125 Å, are 0.237/0.205 meV/GaN below their frozen-cell centerline references. They are not a fully relaxed lower MEP. This is not a whole-path, pressure-relaxed, two-phonon PES or TS certificate. Off-path statics contribute `E+P_ext V` at the same external 45.7 GPa; their residual stress need not vanish. |
| Appendix Fig. A.1, material-path controls | BTO image-count paths and HfO₂ log-strain/linear-cell/CI controls. Source: `figures/vcneb_material_validation.pdf` and the material-validation CSVs. | BTO endpoint rise is 0.08713 eV/BTO, not a barrier. The cited 2.1-kcal/mol BTO comparison is a locally restrained PBEsol DFT distortion cost, **not** a BTO NEB result; that paper's NEB refers to oxygen-vacancy migration. Hf₄O₈ enthalpies are divided by four; literature barriers use different calculator contracts and are not statistical replicates. |
| Table 1, calculator contracts | Per-backend GaN pseudopotential, basis/cutoff, mesh, pressure and force settings. | Do not silently mix different contracts or change the VASP 600-eV path for a prettier comparison. |
| Table 2, acceleration | Matched-start BTO and HfO₂ first crossings of the common 0.10-eV/Å criterion; ABACUS process-launch totals are independently checked against contiguous Slurm job steps. Source: `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md` and `benchmarks/convergence/hf_slurm_abacus_launch_audit_20260928.json`. | BTO 78.5% and HfO₂ 38.2% reductions are case-specific. Slurm timestamps have one-second ties; optimizer records establish the first-crossing grouping. Negative transfers remain archived, and no universal speedup or identical saddle basin is claimed. |

Fig. 8 specifically uses `gan_600eV_atomic_transverse_162_20260929_v2.*`, whose
color is `H(s,q_perp)-H(s,0)`; the archived
`gan_600eV_atomic_dense_surface.*` instead colors absolute `H-H_B4` and is
not the main-paper artwork. The previous 90-point display's negative
interpolant pocket is now resolved by two measured DFT nodes at about
−0.2 meV/GaN; a local frozen-cut lowering is still not a lower relaxed route.

## Archived controls, not additional main figures

The archival 59-point BTO frozen-cubic cut in
`figures/bto_frozen_soft_mode_landscape.*` and its 289-point main-text
Fig. 3 extension are **not** the Fig. 2 restricted sheet. The refined-grid
display-only contour has a 1.61-meV/BTO interior leave-one-out
maximum, and earlier independent coarse-grid checks were worse. It omits stable
atomic motion and cell strain. The 79-point soft-plus-stable frozen cut in
`figures/bto_frozen_path_adapted_79.*` is a different plane again; its
preselected 16-point holdout maximum is 7.19 meV/BTO. Neither is a relaxed
conditional surface. The separate five-point branch-resolved BTO pilot is
recorded in `benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json`;
its local center check does not establish global branch continuity. Do not
merge any of these points into Fig. 2's original nine-node fit.

`figures/gan_600eV_atomic_tube_samples.*` is a 28-point pre-densification
diagnostic without interpolation; its 2.329-meV/GaN maximum linear leave-one-out
error failed the 1-meV gate. The original eight-grid-plus-four-holdout Fig. 5
quadratic model has a 0.04974-meV/GaN axial holdout maximum below its
0.20-meV gate; the final 289-point measured cut remains separate from Fig. 8's path-adapted
chart or a strict TS certificate.

## Reproduction and visual checks

The plotted outputs are PDF/SVG vectors plus PNG previews; 600-dpi TIFFs are
regenerable and repository-ignored. Main multi-panel figures use regular-weight
`(a)`, `(b)` labels, top/right axis spines, readable tick/axis/legend fonts,
aligned panel axes, and semi-transparent white framed legends that do not hide
data. Figure scaling must be checked in the typeset single-document PDF.
The model-colored pixels are never additional DFT calculations; measured
coordinates must stay visible and the caption must name the interpolation.

The dedicated generators are `scripts/plot_cpc_method_and_multibackend_figures.py`
(Figs. 1 and 4), `scripts/plot_bto_qy0_restricted_sheet27.py` (Fig. 2),
`scripts/plot_bto_transverse_soft_landscape.py` (Fig. 3),
`scripts/plot_gan_600eV_local_joint_cut17.py` (Fig. 5),
`scripts/plot_gan_joint_mode_600eV.py` (Fig. 6),
`scripts/plot_gan_gamma_path_manuscript.py` (Fig. 7),
`scripts/plot_gan_600eV_atomic_transverse_q9.py` (Fig. 8), and
`scripts/plot_material_validation_figure.py` (Appendix Fig. A.1). The BTO restricted-sheet
and GaN central-cut generators accept `--output-prefix` and refuse an existing
destination, so
use a fresh prefix for reproduction. Keep source-data and QA JSONs together
with the figures, and never replace manuscript PDF artwork with a stale
root-level export.

The four present BTO/GaN landscape figures have a one-command, DFT-free
rebuild check from committed audit data:
`python -m scripts.rebuild_cpc_landscapes --output-dir tmp/landscape_rebuild`.
Use a fresh output directory. The command requires regenerated source CSVs
to match the archived manuscript tables byte for byte; `--strict-png` also
checks environment-sensitive PNG bytes. The GaN Fig. 5 campaign has now
completed 208 new raw-audited statics and passed its prospective gate; the
rebuild uses the 289-point figure rather than the archived 81-point preview.
On 2026-09-30 all four source CSVs and PNGs regenerated byte-identically;
the new Fig. 5 was also inspected in the compiled manuscript PDF.
