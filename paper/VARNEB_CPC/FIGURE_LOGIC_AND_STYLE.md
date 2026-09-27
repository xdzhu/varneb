# Figure argument and style contract

The five main figures follow one claim sequence: the software keeps control of
the path; its output can be interpreted in a vibrational basis; material paths
preserve their expected topology; the same controller works with distinct
first-principles calculators; and a separate local diagnostic examines
atom--strain coupling at the GaN candidate barrier top. Estimated new-image
evaluation savings appear in Table 3 alongside matched controls and negative
transfers; these totals are not an electronic-launch-by-launch audit.

| Figure | Claim carried by the panels | Source and scientific limit |
| --- | --- | --- |
| 1, architecture | One manager owns the generalized path; worker, calculator, update-policy, cache, and audit concerns have explicit interfaces. | Schematic of the implemented software, not a measured scaling plot. |
| 2, BTO modes | The cubic-endpoint unstable Gamma subspace dominates the atomic T-to-C displacement; enthalpy and strain remain separate coordinates. | Seven-image ABACUS/Phonopy path. Panels (a) mode norms, (b) enthalpy, (c) volume. One deterministic path; no statistical error bars. |
| 3, materials | BTO is monotonic while HfO2 has a resolved interior maximum; image count, interpolation, and CI change distinct aspects of the HfO2 path. | ABACUS source summaries and CSV. Panel (a) uses one BTO formula unit per cell, so its endpoint is **0.08713 eV/f.u.** Panel (b)--(d) divide the Hf4O8 cell by four for energy. The literature barriers use different calculation contracts. |
| 4, GaN backends | Five converged calculators recover the same dominant B4-to-B1 barrier topology, while forward and reverse activation barriers use their own endpoint baselines. | 29 total images at 45.7 GPa, normalized by two GaN units. Panel (a) uses normalized image index; the literature curve is an approximate digitization. Panels (b) and (c) are aligned horizontal-bar comparisons of $H_{\rm peak}-H_{\rm B4}$ and $H_{\rm peak}-H_{\rm B1}$, respectively. Only the forward panel carries the published 0.34-eV/GaN guide; CP2K enters only after the whole path passed audit. |
| 5, GaN local coupling | A local negative direction appears only after atomic and cell-strain coordinates are combined; endpoint optical-mode subspaces geometrically resolve its atomic component. | The two-step joint Hessian is a 1000-eV VASP diagnostic, whereas the endpoint $\Gamma$ bases use the separate 600-eV path protocol. Their projections are geometric overlaps, not energy contributions or a single-protocol barrier. The candidate is near-stationary and connects both endpoint basins locally; neither global path uniqueness nor a finite-temperature transition state follows. Source: `benchmarks/numerical_integrity/gan_ts_endpoint_gamma_bridge_2026-09-27.json`, `figures/gan_joint_mode_coupling_source_data.csv`. |

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
BTO-mode panel with `python scripts/plot_bto_gamma_mode_figure.py --output
paper/VARNEB_CPC/figures/bto_gamma_mode_path`. The archived-summary command
for the material panel is in `docs/material_validation_guide.md`; set its
`--output-dir` to `paper/VARNEB_CPC/figures`. Inspect the final typeset PDF,
not only the standalone PNGs, because figure scaling can make labels too small.
