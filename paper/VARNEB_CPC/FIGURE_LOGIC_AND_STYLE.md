# Figure argument and style contract

The seven main figures follow one claim sequence: the software keeps control
of the path; its output can be interpreted against a vibrational-basis energy
cut; material paths preserve their expected topology; the same controller
works with distinct first-principles calculators; GaN's local atom--strain
coupling, whole-path mode evolution, and directly sampled central enthalpy cut
then expose both the mechanism and its numerical limits. Estimated new-image
evaluation savings appear in Table 2 alongside matched controls and negative
transfers; these totals are not an electronic-launch-by-launch audit.

| Figure | Claim carried by the panels | Source and scientific limit |
| --- | --- | --- |
| 1, architecture | One manager owns the generalized path; worker, calculator, update-policy, cache, and audit concerns have explicit interfaces. | Schematic of the implemented software, not a measured scaling plot. |
| 2, BTO frozen mode cut | A two-soft-direction frozen cubic-cell energy cut and the actual variable-cell path are different objects: the latter projects onto but leaves the plane. | Panel (a): 59 measured ABACUS/PBE points, display-only Clough--Tocher contour and the projected path. Panel (b): seven direct zero-pressure path energies relative to T. Panel (c): off-plane atomic norm; cell strain is omitted too. The contour uses an E−E_C zero, unlike the path's E−E_T zero. The final 59-point interior leave-one-out maximum is 5.14 meV/BTO, while earlier independent coarser-grid holdouts reached 32.57 meV/BTO; no unsampled extremum, conditional PES, or frozen-cut MEP is claimed. Source: `figures/bto_frozen_soft_mode_landscape_source_data.csv` and `_qa.json`. |
| 3, materials | BTO is monotonic while HfO2 has a resolved interior maximum; image count, interpolation, and CI change distinct aspects of the HfO2 path. | ABACUS source summaries and CSV. Panel (a) uses one BTO formula unit per cell, so its endpoint is **0.08713 eV/f.u.** Panel (b)--(d) divide the Hf4O8 cell by four for energy. The literature barriers use different calculation contracts. |
| 4, GaN backends | Five converged calculators recover the same dominant B4-to-B1 barrier topology, while forward and reverse activation barriers use their own endpoint baselines. | 29 total images at 45.7 GPa, normalized by two GaN units. Panel (a) uses normalized image index; the literature curve is an approximate digitization. Panels (b) and (c) are aligned horizontal-bar comparisons of $H_{\rm peak}-H_{\rm B4}$ and $H_{\rm peak}-H_{\rm B1}$, respectively. Only the forward panel carries the published 0.34-eV/GaN guide; CP2K enters only after the whole path passed audit. |
| 5, GaN local coupling | A local negative direction appears only after atomic and cell-strain coordinates are combined; endpoint optical-mode subspaces geometrically resolve its atomic component. | The joint Hessian, endpoint $\Gamma$ bases, and path all use the original 600-eV VASP contract. The 0.02-Å joint step gives one negative eigenvalue, but the 0.0232-eV/Å energy–force mismatch prevents strict TS certification. Optical projections are geometric overlaps, not mode energy contributions. No global path uniqueness or finite-temperature transition-state claim follows. Source: `benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json`, `figures/gan_joint_mode_600eV_source_data.csv`. |
| 6, GaN all-image modes | The full 29-image B4-to-B1 path has compact endpoint-$\Gamma$ atomic projections but substantial independent cell strain. | B4 and B1 optical bases are separate and must not be relabelled as common modes. The leading groups reconstruct this particular mapped atomic path to within 0.00043/0.00139 $\sqrt{\mathrm{amu}}$ Å, while B4-referenced normal strains reach 34.0%, 24.3%, and 17.8%. These are geometric projections, not enthalpy contributions or a globally complete low-dimensional PES. Source: `figures/gan_gamma_path_600eV_source_data.csv` and `_qa.json`. |
| 7, GaN central two-coordinate cut | At a common external pressure of 45.7 GPa, directly sampled $H(s,q_\perp)$ resolves the central barrier region and shows the actual VCNEB centerline within it. | The 90 coordinates comprise 18 path centers and 72 atomic-only transverse statics at their respective fixed cells, all under the original 600-eV VASP contract. Maximum prospective along-path/inner-transverse interpolation errors are 0.166/0.159 meV/GaN under 1-meV gates. This is a bounded frozen central cut, not a whole-path, pointwise pressure-relaxed, conditional-minimum, or TS-certified surface. Source: `figures/gan_600eV_atomic_dense_surface_source_data.csv`, `_qa.json`, and `benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json`. |

The BTO text also reports a separate five-point, branch-resolved conditional
pilot with one preselected holdout. Its source and limitations are recorded in
`benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json`
and `docs/VARNEB_BTO_CONDITIONAL_FIVE_POINT_AUDIT_2026-09-28.md`.
None of those conditional points is silently merged into Fig. 2's 59-point
frozen-cell contour; this pilot does not certify a continuous conditional PES.

The separate BTO path-adapted panel `figures/bto_frozen_path_adapted_79.*`
is a **supplementary candidate, not a replacement for Fig. 2**. It uses a
different plane: $Q_1$ combines the unstable cubic $\Gamma$ modes 0–2 and
$Q_2$ combines stable modes 6–8. Panel (a) shows 79 real fixed-cubic-cell
ABACUS/PBE/100-Ry/10-au-DZP energies and a display-only Clough–Tocher
contour; the seven-image variable-cell path is merely projected. Its T-end
$Q_1$ projection exceeds the sampled grid by
$0.00433\sqrt{\mathrm{amu}}\,$Å, so the plot leaves a white margin rather
than extrapolating the contour. Panel (b) separately shows the direct
variable-cell $H-H_C$ path values, and panel (c) the measured reference-cell
strains. The preselected 16-point interpolation holdout has a maximum
absolute error of 7.19 meV/BTO; this does not certify all unsampled extrema.
The maximum path atomic projection residual is
$0.06028\sqrt{\mathrm{amu}}\,$Å; strain is not part of the frozen plane.
This soft-plus-stable slice must never be relabelled as the main figure's
two-soft-mode plane or its conditional pilot. The source CSVs, source QA,
visual-guide CSV explicitly marked `not_DFT`, and regenerated QA JSON travel
with the supplementary candidate.

An additional GaN central atomic-tube diagnostic, **not** the validated
main-text Fig. 7, is `figures/gan_600eV_atomic_tube_samples.*`. It
places all 28 same-600-eV off-path statics at their measured `(s,q)` values;
no contour or unsampled minimum is drawn. Its maximum linear LOO error
`2.329 meV/GaN` exceeds the predeclared `1.0` gate, so use it only as an
exploratory/supplementary figure until an independent along-path validation
is available. Its source CSV and QA JSON provide per-case hashes.

A distinct **local** same-600-eV joint-coordinate figure is
`figures/gan_600eV_local_joint_cut.*`: panel (a) shows the quadratic
`H(q_u,q_v)` cut only within the measured rectangle and overlays eight raw
grid points plus four axial half-step holdouts; panel (b) compares their
actual DFT enthalpies with model predictions. Its maximum axial holdout
error is `0.04974 meV/GaN`, below the predeclared `0.20` gate. This can
support a local-mode discussion, but must not be captioned as a certified
whole-path surface or transition state. Keep it supplementary unless the
seven-figure manuscript is deliberately restructured and re-typeset.

All plotted panels use a 183-mm figure width, editable SVG/PDF text, and a
600-dpi PNG preview. Panel labels are regular-weight `(a)`, `(b)`, ... at
10.5 pt in the exported figure. Tick labels, axis labels, and legends are
larger than in the initial draft. Plot axes show top and right spines and
ticks; aligned panels use the same grid boundaries. Legends have a white,
partly transparent background and a visible border. Line panels reserve
empty space for the legend or use a shared legend below the grid. No panel
uses an internal title; material, phase, and pressure appear in captions.

Regenerate the architecture and GaN panels from the repository root with
`python scripts/plot_cpc_method_and_multibackend_figures.py`. Regenerate the
BTO frozen-mode panel with the tracked Python generator and its five
hash-matched evidence JSON files:

```text
python scripts/plot_bto_transverse_soft_landscape.py --assembled paper/VARNEB_CPC/evidence/bto_frozen_soft_mode_assembled.json --path-audit paper/VARNEB_CPC/evidence/bto_frozen_soft_mode_path_audit.json --path-report paper/VARNEB_CPC/evidence/bto_frozen_soft_mode_path_report.json --analysis41 paper/VARNEB_CPC/evidence/bto_frozen_soft_mode_analysis41.json --analysis59 paper/VARNEB_CPC/evidence/bto_frozen_soft_mode_analysis59.json --output-prefix <fresh-output-prefix> --allow-exploratory-contours
```

The generator refuses overwrite; use a new output prefix when reproducing.
The earlier pure
mode-norm panel remains archived as `figures/bto_gamma_mode_path.*` and can be
regenerated with `python scripts/plot_bto_gamma_mode_figure.py --output
paper/VARNEB_CPC/figures/bto_gamma_mode_path`. The archived-summary command
for the material panel is in `docs/material_validation_guide.md`; set its
`--output-dir` to `paper/VARNEB_CPC/figures`. Inspect the final typeset PDF,
not only the standalone PNGs, because figure scaling can make labels too small.
