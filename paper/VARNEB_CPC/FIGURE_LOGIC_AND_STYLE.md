# Figure argument and style contract

The four main figures follow one claim sequence: the software keeps control of
the path; its output can be interpreted in a vibrational basis; material paths
preserve their expected topology; and the same controller works with distinct
first-principles calculators. The measured convergence savings appear in
Table 3 alongside their matched controls and negative transfers.

| Figure | Claim carried by the panels | Source and scientific limit |
| --- | --- | --- |
| 1, architecture | One manager owns the generalized path; worker, calculator, update-policy, cache, and audit concerns have explicit interfaces. | Schematic of the implemented software, not a measured scaling plot. |
| 2, BTO modes | The cubic-endpoint unstable Gamma subspace dominates the atomic T-to-C displacement; enthalpy and strain remain separate coordinates. | Seven-image ABACUS/Phonopy path. Panels (a) mode norms, (b) enthalpy, (c) volume. One deterministic path; no statistical error bars. |
| 3, materials | BTO is monotonic while HfO2 has a resolved interior maximum; image count, interpolation, and CI change distinct aspects of the HfO2 path. | ABACUS source summaries and CSV. Panel (a) uses one BTO formula unit per cell, so its endpoint is **0.08713 eV/f.u.** Panel (b)--(d) divide the Hf4O8 cell by four for energy. The literature barriers use different calculation contracts. |
| 4, GaN backends | Five converged calculators recover the same dominant B4-to-B1 barrier topology and similar barrier scale. | 29 total images at 45.7 GPa, normalized by two GaN units. Panel (a) uses normalized image index; the literature curve is an approximate digitization. Panels (b) and (c) show barrier and final maximum generalized force. CP2K enters only after the whole path passed audit. |

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
