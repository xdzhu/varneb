# First complete free-chain update: mechanistic progress, not a final barrier

## Figure argument

Releasing the preserving center exposes transverse structural motion and
lowers its local energy, but neighboring images retain perpendicular
physical force. One free-chain update begins to relax that connection;
it has not produced a converged new barrier.

The four aligned panels show (a) native relative energies, (b) the
transverse y mean-O-minus-Hf displacement proxy, (c) the z proxy associated
with reversal, and (d) ordinary NEB residuals at the seven free interiors.
All data use the same zero-strain partially clamped boundary and physical
ABACUS inputs. The proxies use the common original-T geometry metric,
not electric polarization, mass weighting or phonon occupation.

Blue is the original registered ordinary-converged band. Orange is the
same endpoints and original images with only the released J center
inserted. Green is the complete, native force-evaluated step1 of28949435.
These are three related polygons, not independent statistical replicates
or a comparison of acceleration methods. Every curve has its own actual
normalized generalized arc. Straight connectors are not a fitted smooth
MEP. The seed/step1 curves overlap closely because this is just one update;
no visual offsets, smoothing or synthetic points separate them.

| Polygon | Sampled peak (meV/HfO2) | NEB fmax (eV/A) | Limiting image |
| --- | ---: | ---: | ---: |
| Original | 96.828183 | 0.098817 | 4 |
| Released-center seed | 85.499889 | 0.125611 | 3 |
| Free-chain step1 | 85.269291 | 0.123740 | 3 |

For the limiting vector at step1, perpendicular and spring projections
are+0.125157 and-0.001417eV/A, summing to the residual0.123740. The
perpendicular-only maximum is0.125219, whereas the spring-only maximum
is0.007727. These are instantaneous generalized-force components, not
energies or a prediction of what changing k would achieve. k stays0.2.
The large force is not simply a spring-spacing artifact. Source data
retain the signed components and their verified vector-sum defect.

`source_data.csv` has27 actual image rows; endpoint NEB residuals are
blank because those two images are fixed, not imputed zero-force samples.
The selected new frame has seven fresh SCFs,6.557367native-proxy core-hours.
All seed/local/preparation and other completed-frame costs are additional.
The archive's9 native E/F/stress records were independently reparsed
locally. Undistributed pseudopotential/orbital bytes are not claimed to
have been rehashed locally. No new DFT or scheduler mutation is invoked.

## Reproduce

```text
python scripts/plot_hfo2_preserving_connection.py --step 1 --output benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/figures/step_0001
python -m pytest -q tests/test_hfo2_preserving_connection_snapshot.py tests/test_hfo2_preserving_connection_figure.py
```

PDF/SVG retain editable text; PNG is the300dpi preview. Normal13pt(a)--(d)
labels, larger fonts, top/right spines, aligned axes and white85%-opaque
framed legends follow the author's existing style. The preview was
visually checked for clipping, overlaps and misleading curve separation.
This progress figure stays in the case archive rather than replacing a
main-paper final-result figure before the connection test finishes.
