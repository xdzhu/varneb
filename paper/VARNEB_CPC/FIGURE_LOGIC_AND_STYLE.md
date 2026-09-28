# VARNEB CPC figure argument and evidence contract

This file describes the compiled two-column manuscript as of 2026-09-28:
**five main figures, four supplementary figures, and two main tables**. Figure
numbers below follow `varneb_CPC.tex`; S1–S4 follow
`varneb_CPC_supplement.tex`. Measured DFT markers, fitted contours, projected
paths, and external literature curves must remain visually and verbally distinct.

| Item | Claim and source | Boundary that the caption and text must retain |
| --- | --- | --- |
| Fig. 1, architecture | One manager owns the generalized path; workers, calculators, update policies, cache, and audits have separate interfaces. Source: `figures/varneb_architecture.pdf`. | An implementation diagram, not a measured scaling benchmark. |
| Fig. 2, BTO restricted sheet | ABACUS/PBE/100-Ry/10-au-DZP BTO at zero pressure. Nine raw-audited model-training nodes and two *prospective* DFT checks define a local `Q_y=0` variable-cell sheet; the seven-image T→C path is projected onto it. Sources: `figures/bto_qy0_restricted_even_mode_sheet_points_source_data.csv`, `_contour_source_data.csv`, and `_qa.json`. | Contour is the sixth-order even-mode model only inside the measured-coordinate hull. The two blind energy errors are −0.150 and −0.977 meV/BTO against a ±2-meV gate. Four older checks and the separately optimized T endpoint are retrospective, not three more blind points. Stable atomic modes and strain are released but not the plotted axes; the omitted cubic soft mode is fixed at `Q_y=0`. Neither an unrestricted conditional PES, finite-temperature FES, nor T→C activation barrier. |
| Fig. 3, GaN backends | Five converged 29-total-image B4→B1 paths at 45.7 GPa, with forward/reverse barriers referenced to their respective endpoints. Source: `figures/gan_multibackend_validation.pdf` and its committed source tables. The CP2K all-interior raw-energy/stress versus exact-cache report is `benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json`. | Energies are enthalpies per GaN. Literature curve is approximate digitization, not matched raw data; only the forward panel uses the published 0.34-eV/GaN guide. Calculator contracts differ, so cross-code agreement is topological and approximate, not identity. CP2K raw atomic forces and final run-end markers remain unavailable, despite 27/27 energy/stress reconciliation. |
| Fig. 4, GaN local joint cut | Original VASP/PBE 600-eV joint atom–strain coordinates show a local negative and a positive curvature near the highest image. The final display interpolates 81 measured enthalpies on a complete 9×9 grid. Sources: `figures/gan_600eV_local_joint_dft81_v2_20260928_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_ts_2d_9x9_refinement_20260928.json`. | The prior 25-point interpolant predicted 56 new nodes within 0.0135 meV/GaN; the original quadratic model's ±0.2-meV/GaN gate applied only to four predeclared axial holdouts. This frozen local cut is neither an orthogonally relaxed or whole-path surface nor a strict TS certificate; the 0.0232-eV/Å energy–force mismatch remains. |
| Fig. 5, GaN central 2D cut | The 45.7-GPa central `H(s,q_perp)-H(s,0)` contour resolves the small transverse scale separately from the full path barrier in panel b. It uses the same 18 measured path centers and 72 atomic-only off-path statics under the original VASP/PBE 600-eV contract. Sources: `figures/gan_600eV_atomic_dense_surface_source_data.csv`, `gan_600eV_atomic_transverse_landscape_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json`. | The interpolant is linear in path coordinate and quadratic in transverse displacement only within the 90-point sampled rectangle. Prospective along-path/inner-transverse maximum errors are 0.166/0.159 meV/GaN, under 1-meV gates. This is a frozen *central* cut, not a whole-path, pressure-relaxed, two-phonon PES or TS certificate. The off-path statics contribute `E+P_ext V` at the same external 45.7 GPa; their residual stress need not vanish. |
| Supplementary Fig. S1, material-path controls | BTO image-count paths and HfO₂ log-strain/linear-cell/CI controls. Source: `figures/vcneb_material_validation.pdf` and the material-validation CSVs. | BTO endpoint rise is 0.08713 eV/BTO, not a barrier. The cited 2.1-kcal/mol BTO comparison is a locally restrained PBEsol DFT distortion cost, **not** a BTO NEB result; that paper's NEB refers to oxygen-vacancy migration. Hf₄O₈ enthalpies are divided by four; literature barriers use different calculator contracts and are not statistical replicates. |
| Supplementary Fig. S2, whole-path GaN modes | Two independent endpoint Γ optical bases, B4-referenced strain, and reconstruction residuals across the mapped 29-image path. Sources: `figures/gan_gamma_path_600eV_source_data.csv` and `_qa.json`. | Five vectors in three leading optical groups capture >99.99% of this path's squared optical displacement; the B4/B1 bases are not mode-to-mode identified. Maximum residuals are 0.00043/0.00139 `sqrt(amu) Å`, while normal strains reach 34.0%, 24.3%, and 17.8%. This is geometry, not a decomposition of enthalpy. |
| Supplementary Fig. S3, GaN Hessian and Γ projections | Frozen atomic, frozen strain, joint and strain-relaxed curvatures accompany the local Fig. 4 cut; endpoint optical-group projections describe the joint negative direction. Sources: `figures/gan_joint_mode_600eV_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json`. | The coordinate metric controls numerical curvature magnitudes. Endpoint optical overlaps are geometric, not mode-resolved barrier energies; the joint candidate is not a certified TS. |
| Supplementary Fig. S4, BTO frozen two-soft-mode grid | The original 59 static samples were extended by 22 new raw-audited ABACUS/PBE/100-Ry/10-au-DZP points to a complete 9×9 grid. Sources: `figures/bto_frozen_soft_mode_81_20260928_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/bto_transverse_soft_frozen81_20260928.json`. | The prior 59-point interpolant predicts the new nodes within 0.441 meV/BTO, but full-grid interior LOO maximum is 5.14 meV/BTO and the display interpolant overshoots a measured minimum by 1.78 meV/BTO. The projected VCNEB path leaves the frozen plane; this is not the main Fig. 2 restricted sheet, a relaxed PES, or a barrier. |
| Table 1, calculator contracts | Per-backend GaN pseudopotential, basis/cutoff, mesh, pressure and force settings. | Do not silently mix different contracts or change the VASP 600-eV path for a prettier comparison. |
| Table 2, acceleration | Matched-start BTO and HfO₂ first crossings of the common 0.10-eV/Å criterion; ABACUS process-launch totals are independently checked against contiguous Slurm job steps. Source: `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md` and `benchmarks/convergence/hf_slurm_abacus_launch_audit_20260928.json`. | BTO 78.5% and HfO₂ 38.2% reductions are case-specific. Slurm timestamps have one-second ties; optimizer records establish the first-crossing grouping. Negative transfers remain archived, and no universal speedup or identical saddle basin is claimed. |

Fig. 5 specifically uses `gan_600eV_atomic_transverse_landscape.*`, whose
color is `H(s,q_perp)-H(s,0)`; the archived
`gan_600eV_atomic_dense_surface.*` instead colors absolute `H-H_B4` and is
not the main-paper artwork. All sampled off-path DFT excess enthalpies in Fig. 5
are nonnegative. The bounded quadratic interpolant alone dips to
−0.262 meV/GaN between nodes; this is not evidence of a lower DFT route.

## Archived controls, not additional main figures

The archival 59-point BTO frozen-cubic cut in
`figures/bto_frozen_soft_mode_landscape.*` and its 81-point Supplementary
Fig. S4 extension are **not** the Fig. 2 restricted sheet. The full-grid
display-only contour has a 5.14-meV/BTO interior leave-one-out
maximum, and earlier independent coarse-grid checks were worse. It omits stable
atomic motion and cell strain. The 79-point soft-plus-stable frozen cut in
`figures/bto_frozen_path_adapted_79.*` is a different plane again; its
preselected 16-point holdout maximum is 7.19 meV/BTO. Neither is a relaxed
conditional surface. The separate five-point branch-resolved BTO pilot is
recorded in `benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json`;
its local center check does not establish global branch continuity. Do not
merge any of these points into Fig. 2's nine-node fit.

`figures/gan_600eV_atomic_tube_samples.*` is a 28-point pre-densification
diagnostic without interpolation; its 2.329-meV/GaN maximum linear leave-one-out
error failed the 1-meV gate. The original eight-grid-plus-four-holdout Fig. 4
quadratic model has a 0.04974-meV/GaN axial holdout maximum below its
0.20-meV gate; the final 81-point measured cut remains separate from Fig. 5's path-adapted
chart or a strict TS certificate.

## Reproduction and visual checks

The plotted outputs are PDF/SVG vectors plus PNG previews; 600-dpi TIFFs are
regenerable and repository-ignored. Main multi-panel figures use regular-weight
`(a)`, `(b)` labels, top/right axis spines, readable tick/axis/legend fonts,
aligned panel axes, and semi-transparent white framed legends that do not hide
data. Figure scaling must be checked in the typeset main and supplement PDFs.
The model-colored pixels are never additional DFT calculations; measured
coordinates must stay visible and the caption must name the interpolation.

The dedicated generators are `scripts/plot_cpc_method_and_multibackend_figures.py`
(Figs. 1 and 3), `scripts/plot_bto_qy0_restricted_sheet.py` (Fig. 2),
`scripts/plot_gan_600eV_local_joint_cut.py` (Fig. 4),
`scripts/plot_gan_600eV_atomic_transverse_landscape.py` (Fig. 5; same sampled source as the absolute-contour audit),
`scripts/plot_material_validation_figure.py` (S1), and
`scripts/plot_gan_gamma_path_manuscript.py` (S2), and
`scripts/plot_gan_joint_mode_600eV.py` (S3). The BTO restricted-sheet
generator verifies evidence hashes but writes its named figure/source files
in place; do not run it over a dirty figure checkout. The GaN central-cut
generator accepts `--output-prefix` and refuses an existing destination, so
use a fresh prefix for reproduction. Keep source-data and QA JSONs together
with the figures, and never replace manuscript PDF artwork with a stale
root-level export.
