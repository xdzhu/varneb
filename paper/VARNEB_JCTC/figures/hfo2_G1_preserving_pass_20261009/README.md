# Figure contract: an ordinary pass is not a stationary bottleneck

Conclusion: the preserving flip has reached ordinary NEB 0.10 eV/A and has
two sampled maxima flanking a lower Pbcn centre, but physical tangential
forces at those maxima still preclude stationary-TS certification.

This dated quantitative grid uses the existing Python/matplotlib workflow.
It does not replace the earlier 6/20/14/10 observation figure. Four actual
complete source snapshots, 6/39/69/32, are replayed against frozen descriptors
before drawing; changed steps/pass labels or claimed G1/TS gates are rejected.

## Caption

Complete same-PO+ observations. (a) Low-energy discrete profiles for PO-to-T,
PO-to-M and the preserving flip. (b) The still-unconverged reversing profile
on its own energy scale, with no clipping or hidden high-energy values.
(c) Ordinary moving-image NEB residuals; the dashed line is 0.10 eV/A.
(d) Preserving physical tangential-force magnitude versus atomic/cell NEB
residuals in the same original generalized-coordinate metric. The ordinary
threshold is not an independently established physical-tangent TS tolerance.
Fixed endpoints have no NEB residual. Connections join calculated samples,
not a validated smooth MEP, a mode energy decomposition or a conditional PES.

All energies are relative to identical PO+, -9783.249675811956 eV/Hf4O8,
divided by four formula units; P=0, E=0. Normalized source arcs preserve their
original metrics and lifts; this normalization is not a common physical length.
The Pbcn annotation identifies the central image4, **not** either peak.
Reference geometric labels are not assigned phonons or electronic polarization.
No transverse-stability, full-index or sampling-error certificate is supplied.

## Export, statistics and visual checks

183x177 mm; editable SVG/PDF text, PNG preview and 600-dpi LZW TIFF. Normal
13-pt (a)--(d), no panel titles, visible top/right spines and inward ticks.
Axes/panel-label rows and columns align. Shared framed white legend and
separate force legend sit outside curves with alpha0.85; labels/ticks are
enlarged to follow the user's explicit style. Source CSV includes all37
records including cached/reused endpoints; they are not independent replicates.
No inferred statistical error bars, tests or p-values are shown. No raster
retouching, smoothing, contrast/crop manipulation or external-paper image reuse.

The generated QA records numeric/layout checks and raw export hashes. Actual
pixel and compiled-PDF review are separately recorded in the delivery receipt;
the pending review flag in the original generated QA is not retroactively edited.

```sh
python -m scripts.plot_hfo2_G1_preserving_pass \
  --specification benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/network_specification.json \
  --material-report benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/network_progress_local.json \
  --chain-report benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/analysis_local.json \
  --output /new/path/figure
python -m pytest -q tests/test_hfo2_preserving_pass_figure.py
```

The entry refuses an existing output before reading input. Original observations,
symmetry sensitivity, calculator/SCF checks and source identities are in the
[material case](../../../../benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/README.md).
