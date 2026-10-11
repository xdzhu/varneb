# Clamped PO-to-M: figure contract

Conclusion: along the first ordinary-converged common-substrate escape chain,
the discrete energy maximum, reference-space atomic coverage and allowed cell
shear describe different aspects of the transformation; a selected atomic
projection alone does not reconstruct the joint path.

Quantitative grid, Python/matplotlib, following the existing repository plotting
workflow. Panel (a) is the discrete P=0 energy path; (b) is squared atomic
displacement coverage in the common original-T mass metric; (c) is finite Green
strain in the original-T axes. All share the original normalized generalized
arc. The dashed guide locates sampled image 3, not a certified TS. Nine native
images, including two cached endpoints, are not independent statistical repeats.

Sources: frozen E068 terminal observation and E071 corrected-frame numeric
analysis. Their identity is checked before plotting. This script does not repeat
the native SCF audit or redistribute author mode arrays. Four Cmma registrations
are all retained; their narrow envelope is an orientation range, not an error
bar. Rank three versus rank four is descriptive, not a fair model-accuracy test.
No energy partition, local phonon, conditional PES or forecast is inferred.

Original-T x/z form the fixed plane; original-T y is the open vector. Green
Eyy, Exy and Eyz are dimensionless tensor entries, multiplied by 100 for display,
not engineering shear or imposed strain increments. The imposed training strain
is zero, even though the open vector can shear. Atomic coverage excludes these
cell coordinates.

183 x 207 mm, editable SVG/PDF, 300 dpi PNG and 600 dpi LZW TIFF; source CSV and
machine-readable QA. Normal-weight larger (a) labels; aligned axes; no panel
titles; top/right spines and ticks; framed translucent white legends outside the
data panels. Straight connectors only; no interpolation or smoothing. Actual
visual QA is required after generation. No new DFT, selection or holdout access.

Run from the repository root with a *fresh* output directory (the archived bundle
is not overwritten):

```text
python -m scripts.plot_hfo2_clamped_M --output NEW_DIRECTORY
```

Actual PNG review finds no clipping or legend/curve overlap, and the SVG/PDF
contain editable text. The complete Git archive reproduces the same nine CSV
rows using the same installed Python/library environment; the E071 JSON content
pin deliberately tolerates checkout line endings while preserving all values.
Raster/vector byte identity across different render runs is not claimed.
Local and genuinely archive-rooted tests pass 28 each with direct exit0. A
preliminary command labelled "clean" was run from the working tree and is not
counted as independent archive verification. See validation.json for the
actual tested tree, manuscript compilation and scheduler-observation scope.
