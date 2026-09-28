# GaN local landscape promotion and rendered-manuscript QA — 2026-09-28

## Figure contract

- Core conclusion: The highest-image GaN candidate has a resolved *local*
  atom--strain-coupled negative direction and a positive joint-transverse
  direction under the original 600-eV VASP contract at 45.7 GPa.
- Archetype: quantitative two-panel grid, with the joint-coordinate contour
  as the hero and independently measured model checks as validation.
- Backend: the existing Python/matplotlib figure generator; no new rendering
  code, DFT evaluation, electronic-structure setting, or source data were
  introduced by the manuscript promotion.
- Output: the existing vector PDF in main Fig. 4; source CSV, QA JSON, SVG,
  and PNG remain paired with it in `figures/`.
- Reviewer risk: the frozen quadratic rectangle does not imply a whole-path
  two-mode PES, orthogonally relaxed surface, or certified variable-cell TS.
  The 0.0232-eV/Å energy--force discrepancy remains explicit in the text.

## Evidence-to-panel map

| Panel | Evidence | Limit |
| --- | --- | --- |
| Main Fig. 4a | Eight measured grid statics, four preselected axial half-step statics, and the highest-image center; original 600-eV/PBE/Ga_d+N VASP contract at 45.7 GPa. | Quadratic interpolation only inside the sampled `q_u/q_v` rectangle; no orthogonal relaxation. |
| Main Fig. 4b | DFT--model enthalpy comparison; largest axial holdout error 0.04974 meV/GaN against a predeclared 0.20-meV/GaN gate. | Axial checks do not validate every off-axis interior point. |
| Supplementary Fig. S3a | Original atom, strain, joint, and strain-relaxed Hessian-block chart. | Eigenvalue magnitudes depend on the stated joint coordinate metric. |
| Supplementary Fig. S3b | Endpoint-Γ optical-group projections of the joint negative direction. | Geometric overlap, not a mode-resolved barrier-energy partition. |

The scientific source audit is `benchmarks/numerical_integrity/gan_600eV_ts_2d_pilot_20260928.json`;
the corresponding figure records are
`figures/gan_600eV_local_joint_cut_source_data.csv`,
`figures/gan_600eV_local_joint_cut_qa.json`, and
`figures/gan_joint_mode_600eV_qa.json`. Main Fig. 5 remains the separately
validated 90-point *central* path-adapted frozen cut; its coordinates and
energy zero are not merged with Fig. 4.

The attractive 59-point BTO frozen-cubic two-soft-mode contour remains an
archived exploratory visual, not a replacement for main Fig. 2. Its projected
T-to-C path leaves the frozen plane through stable-mode motion and cell
strain; independent coarse-grid errors were as large as 32.57 meV/BTO.
Main Fig. 2 instead uses the measured-hull-limited, `Q_y=0` restricted
variable-cell sheet and two prospective sub-2-meV/BTO checks. The main
BTO caption and text explicitly distinguish a coordinate projection from
the actual variable-cell path.

## Compile and visual checks

The TeX Live/latexmk build returned success for both documents, with no
undefined citations/references or LaTeX errors. The main paper remains 10
pages; the supplement is three pages. The only overfull warning is the
pre-existing 1.9-pt box on main page 1. Main page 8 and supplement page 3
were rendered from the resulting PDFs and inspected at 1.8× scale: panel
labels, axes, markers, model colorbar, and captions remain readable and
unclipped. The additional supplement page does not alter the main-paper
page count. Author metadata and final journal-specific export checks remain
open submission gates.
