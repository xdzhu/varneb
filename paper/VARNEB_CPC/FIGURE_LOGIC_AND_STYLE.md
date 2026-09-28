# VARNEB CPC figure argument and evidence contract

This file describes the compiled two-column manuscript as of 2026-09-28:
**five main figures, two supplementary figures, and two main tables**. Figure
numbers below follow `varneb_CPC.tex`; S1–S2 follow
`varneb_CPC_supplement.tex`. Measured DFT markers, fitted contours, projected
paths, and external literature curves must remain visually and verbally distinct.

| Item | Claim and source | Boundary that the caption and text must retain |
| --- | --- | --- |
| Fig. 1, architecture | One manager owns the generalized path; workers, calculators, update policies, cache, and audits have separate interfaces. Source: `figures/varneb_architecture.pdf`. | An implementation diagram, not a measured scaling benchmark. |
| Fig. 2, BTO restricted sheet | ABACUS/PBE/100-Ry/10-au-DZP BTO at zero pressure. Nine raw-audited model-training nodes and two *prospective* DFT checks define a local `Q_y=0` variable-cell sheet; the seven-image T→C path is projected onto it. Sources: `figures/bto_qy0_restricted_even_mode_sheet_points_source_data.csv`, `_contour_source_data.csv`, and `_qa.json`. | Contour is the sixth-order even-mode model only inside the measured-coordinate hull. The two blind energy errors are −0.150 and −0.977 meV/BTO against a ±2-meV gate. Four older checks and the separately optimized T endpoint are retrospective, not three more blind points. Stable atomic modes and strain are released but not the plotted axes; the omitted cubic soft mode is fixed at `Q_y=0`. Neither an unrestricted conditional PES, finite-temperature FES, nor T→C activation barrier. |
| Fig. 3, GaN backends | Five converged 29-total-image B4→B1 paths at 45.7 GPa, with forward/reverse barriers referenced to their respective endpoints. Source: `figures/gan_multibackend_validation.pdf` and its committed source tables. | Energies are enthalpies per GaN. Literature curve is approximate digitization, not matched raw data; only the forward panel uses the published 0.34-eV/GaN guide. Calculator contracts differ, so cross-code agreement is topological and approximate, not identity. |
| Fig. 4, GaN local coupling | Original VASP/PBE 600-eV path, endpoint Γ bases, and atom–strain finite differences show a local coupled negative direction. Sources: `figures/gan_joint_mode_600eV_source_data.csv`, `benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json`. | The 0.0232-eV/Å energy–force mismatch prevents a strict full-variable-cell TS certificate. Endpoint optical subspaces are geometric projections, not mode-resolved barrier energies. |
| Fig. 5, GaN central 2D cut | The 45.7-GPa central `H(s,q_perp)` chart has 18 measured path centers and 72 atomic-only off-path statics, all under the original VASP/PBE 600-eV contract. Sources: `figures/gan_600eV_atomic_dense_surface_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json`. | The interpolant is linear in path coordinate and quadratic in transverse displacement only within the 90-point sampled rectangle. Prospective along-path/inner-transverse maximum errors are 0.166/0.159 meV/GaN, under 1-meV gates. This is a frozen *central* cut, not a whole-path, pressure-relaxed, two-phonon PES or TS certificate. The off-path statics contribute `E+P_ext V` at the same external 45.7 GPa; their residual stress need not vanish. |
| Supplementary Fig. S1, material-path controls | BTO image-count paths and HfO₂ log-strain/linear-cell/CI controls. Source: `figures/vcneb_material_validation.pdf` and the material-validation CSVs. | BTO endpoint rise is 0.08713 eV/BTO, not a barrier. Hf₄O₈ enthalpies are divided by four; literature barriers use different calculator contracts and are not statistical replicates. |
| Supplementary Fig. S2, whole-path GaN modes | Two independent endpoint Γ optical bases, B4-referenced strain, and reconstruction residuals across the mapped 29-image path. Sources: `figures/gan_gamma_path_600eV_source_data.csv` and `_qa.json`. | Five vectors in three leading optical groups capture >99.99% of this path's squared optical displacement; the B4/B1 bases are not mode-to-mode identified. Maximum residuals are 0.00043/0.00139 `sqrt(amu) Å`, while normal strains reach 34.0%, 24.3%, and 17.8%. This is geometry, not a decomposition of enthalpy. |
| Table 1, calculator contracts | Per-backend GaN pseudopotential, basis/cutoff, mesh, pressure and force settings. | Do not silently mix different contracts or change the VASP 600-eV path for a prettier comparison. |
| Table 2, acceleration | Matched-start BTO and HfO₂ first crossings of the common 0.10-eV/Å criterion; ABACUS process-launch totals are independently checked against contiguous Slurm job steps. Source: `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md` and `benchmarks/convergence/hf_slurm_abacus_launch_audit_20260928.json`. | BTO 78.5% and HfO₂ 38.2% reductions are case-specific. Slurm timestamps have one-second ties; optimizer records establish the first-crossing grouping. Negative transfers remain archived, and no universal speedup or identical saddle basin is claimed. |

## Archived controls, not additional main figures

The 59-point BTO frozen-cubic two-soft-mode cut in
`figures/bto_frozen_soft_mode_landscape.*` is **not** the Fig. 2 restricted
sheet. Its display-only contour had a 5.14-meV/BTO final interior leave-one-out
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
error failed the 1-meV gate. The eight-grid-plus-four-holdout local quadratic
`figures/gan_600eV_local_joint_cut.*` has a 0.04974-meV/GaN axial holdout
maximum below its 0.20-meV gate, but is a separate local model, not Fig. 5's
path-adapted chart or a strict TS certificate.

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
`scripts/plot_gan_joint_mode_600eV.py` (Fig. 4),
`scripts/plot_gan_600eV_atomic_dense_surface.py` (Fig. 5),
`scripts/plot_material_validation_figure.py` (S1), and
`scripts/plot_gan_gamma_path_manuscript.py` (S2). The BTO restricted-sheet
generator verifies evidence hashes but writes its named figure/source files
in place; do not run it over a dirty figure checkout. The GaN central-cut
generator accepts `--output-prefix` and refuses an existing destination, so
use a fresh prefix for reproduction. Keep source-data and QA JSONs together
with the figures, and never replace manuscript PDF artwork with a stale
root-level export.
