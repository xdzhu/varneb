# Reconnected ordinary preserving band at zero strain

The figure compares three actual polygons with the same ordered atoms,
fixed endpoint caches, partial clamp and electronic inputs. Panel (a)
shows the sampled energies, (b) the transverse mean O-minus-Hf structural
displacement, (c) its z component, and (d) the seven interior NEB residuals.
The structural proxies are not electronic polarization. Every series uses
its own normalized generalized arc; lines are straight connectors, not a
continuous MEP fit. Endpoint residuals are undefined and omitted, not zero.

Original job28949435 completed0:0 after9updates/63fresh interior SCFs.
The independently reconstructed final residual is0.0940120033eV/A and
the sampled maximum is83.1662809401meV/HfO2, compared with96.8281830028
for the original band. The central transverse displacement remains
0.1210744444Angstrom. All seven interiors were free, without the four
local-intervention constraints; the same endpoint states were retained.

This establishes a lower connected **sampled ordinary-NEB candidate**,
not a certified continuous saddle or global MEP. Its1.250meV/HfO2
separation from the original zero-strain T candidate is not resolved as
a robust ordering by the current sampling and force target. The original
registered training labels and independent forecast scores are unchanged.
One seed was tested, not statistical replicates or an acceleration benchmark.

The complete519-file native archive retains all63new calls, exact
INPUT/KPT/STRU and native logs. All energy/force/stress results were
reparsed independently locally; full pseudo/orbital bytes were checked
on HF, not redistributed or rehashed locally. Runtime costs are55.718818
native-call proxy and55.857778allocation core-hours; earlier preparation,
local release and original-band costs are additional.

Python/matplotlib generated editable PDF/SVG, PNG,27-row source CSV and
QA metadata. No smoothing, synthetic observations, curve offsets or
statistical error bars were applied. Normal13pt panel labels, visible
top/right axes, larger fonts and framed semitransparent white legends
follow the author style. The PNG was visually reviewed at final layout.

```text
python scripts/plot_hfo2_preserving_connection.py --terminal --step 9 --output benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/figures/terminal_step_0009
```
