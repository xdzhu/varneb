# Mechanism-focused main figures

Python/matplotlib, using the existing calculation and geometry-analysis workflow.
No new DFT, refitted modes, atom remapping, physical input changes or unseen-strain
labels are used. The existing eight-band result and ordered endpoint/structure
records are inputs, not a new whole-campaign audit.

## Visual argument

1. **Competition:** a decrease of the high reversing candidate does not lower
   the easiest examined switching peak. The observed channel responses are
   different from their final-well-only predictions. Two quantitative panels:
   paired sampled peaks and actual versus endpoint-only response.
2. **Path geometry:** the same sign reversal of the polar displacement can
   accompany preservation or reversal of an oxygen shuffle. Compare the two
   flips and T escape using signed geometric coordinates and open-cell strain,
   with both substrate conditions. These are structural descriptions, not
   phonon energies, Berry polarization or an isolated causal energy partition.
3. **Matched relaxation:** at the fixed retained coordinates of the two local
   cohorts, atomic relaxation and additional joint relaxation have different
   measured energy scales. Retain all five coordinates in each cohort and
   do not reinterpret their differences as whole-path barrier corrections.

Figures use quantitative grids with consistent condition/treatment colors,
183 mm width, editable PDF/SVG text and PNG previews. Labels `(a)` etc. have
normal weight; top/right frames remain; axes and panel labels align. Legends
have translucent white framed backgrounds and are outside data panels.
No individual panel titles, invented uncertainty bars, smoothing or missing
data imputation. The molecular coordinates use the fixed original T reference
metric and an initial-image periodic gauge held fixed for each whole chain.

Source CSV/JSON files and a single replay script accompany the plots. Samples
are deterministic images or matched calculation points, not independent
statistical replicates. Main-text interpretation is bounded to the registered
candidate paths and prescribed mechanical ensemble.

## Interpretation and visual review, 2026-10-18

The signed Qx pattern distinguishes preserving from reversing even though
both reverse the polar-displacement proxy. T escape replaces the initial
Qx pattern by Qz. The remaining atomic displacement is exported explicitly:
near-zero triplet amplitudes do not identify a cubic structure. The response
plot also removes the common initial-well shift, preventing that offset
from being mistaken for the cause of opposite switching responses.

Local treatment bars demonstrate atomic relaxation at each fixed seed cell;
the orange difference includes further atomic readjustment and is not a
unique elastic energy. The Q1/Q2 axes in this treatment are joint coordinates,
not the Qx/Qz geometric projection. All energies are per HfO2 formula unit.
The different energy ranges for the reversing column and the two treatment
panels are declared in the manuscript captions.

All three exported PNGs were visually inspected. The LuaLaTeX reading PDF
compiled successfully to 46 pages, including the retained appendices.
All-page contact sheets were inspected for missing/overflowing content;
main-figure pages 10, 14 and 15 were inspected at reading resolution.
No clipped labels, covered data, detached figure captions or overlapping
panel axes were found. The compiler reported no overfull/underfull boxes or
missing characters. The 130 scoped geometry, plotting and related manuscript
tests passed; this is not a full repository test run or a new physical audit.

The original nine figure assets, all 213 table lines and external citations
are retained. Main Results plus Discussion decreased from 4478 to 2979 words
while adding the mechanism interpretation. No complete article or new-method
validation is claimed: the separate +0.5% path labels are not included.
