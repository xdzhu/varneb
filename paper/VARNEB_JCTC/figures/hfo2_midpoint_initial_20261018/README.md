# Same-target midpoint seed comparison: initial snapshots, not final paths

## Scientific question and evidence contract

At the identical +0.5% substrate and exact same endpoint caches, do the
zero-strain and +1%-strain initial polygons differ predominantly around
preserving reversal, T escape, or both? Energy differences are shown
alongside oxygen rearrangement and allowed-cell response, not assigned
uniquely to one of those coordinates.

The quantitative six-panel grid uses the established Python/matplotlib paper
workflow at183mm width. (a,b) show energies, (c,d) signed geometric Qx/Qz,
and(e,f) original-T open-y stretches. The left column is preserving reversal,
right is PO-to-T. Blue is the completed E268 ordinary chain; orange is
only the complete E269 step-zero evaluation. The legend explicitly distinguishes
those unequal optimization stages. All curves are at the same target strain;
0% and1% identify their initializing polygons, not different plotted substrates.

## Measured stage result

| Channel | E268 ordinary peak | E269 initial peak | Initial minus E268 |
| --- | ---: | ---: | ---: |
| preserving | 98.530877 | 138.889103 | 40.358226 |
| PO-to-T | 90.418178 | 90.503363 | 0.085185 |

Energies are meV/HfO2, four formula units per twelve-atom cell. The orange
chains do not pass0.10eV/A: residuals are0.135999 and0.151237eV/A. Both original
jobs continue optimizing. These are not final start sensitivities, an error
bar on a barrier, new independent predictions or new saddle certificates.
The original holdout labels and scores remain unchanged.

Preserving shows a concurrent structural difference: peak Qx is -0.798917
versus -1.039388A and open-y stretch -0.416145 versus+0.654708%. Both remain
in the preserving sign sector. This does not isolate an atomic or cell cause:
both structures and cell vectors differ. T profiles are energetically close
already, although the orange T chain is not ordinary-converged.

## Integrity and visual checks

Nine actual images per curve, including the same cached endpoints; deterministic
configurations, not independent statistical repeats. Every line uses its own
normalized generalized arc. Straight connectors are not a smooth optimized MEP.
No image remapping, artificial error bars, nonlinear fitting or curve smoothing.
Qx/Qz are geometric oxygen patterns, not phonon occupations or assigned irreps.

PDF/SVG text remains editable; PNG preview is300dpi. Normal-weight13pt(a) labels,
top/right frames, aligned axes and a shared translucent-white framed legend follow
the author's style. The final PNG was visually checked. The first export put a
secondary legend close to(b); the final version uses one shared two-row legend,
clear of all labels and data. There are no clipped labels or annotations on curves.

The38 scoped tests cover source geometry, native initial data replay, units,
the refusal to score a running/failed job as terminal, arithmetic, and figure
sample/label counts. They do not certify physical convergence.

Source data are comparison.json, outcomes.csv and profiles.csv. Native E269
step-zero observations are in the dated case's initial_complete_chains folder.
Regenerate from the repository root with:

```sh
python benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/compare_paths.py --exported benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/initial_complete_chains --output NEW_FIGURE_DIRECTORY --initial --draw
```

This is a research-stage figure, not a replacement for the manuscript's final
ordinary-path comparison. Completed E269 results must be processed separately.
