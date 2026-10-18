# Independent physical-strain test and mechanism figures

The two +0.5% paths complete the preregistered sampled-peak test. Preserving
switching and PO-to-T have peaks 98.530877 and 90.418178 meV/HfO2, with
NEB residuals 0.079773 and 0.080825 eV/A. No physical input was changed.

The selected-pair contrast T minus preserving changes from -12.411559 at
zero strain to -8.112699 at +0.5%, then -28.660781 at +1%. The intermediate
condition improves this contrast but does not reverse its sign. Only two
channels were tested there, so it does not establish the four-edge network.

B0 midpoint interpolation overpredicts preserving by12.595430 meV/HfO2 but
T by only0.171959. Endpoint-only B1 predictions are closer for preserving
and further from T. All138 frozen records remain in analysis.json and
forecast_errors.csv:6numeric values and132 unavailable transfers. This
does not demonstrate a physical prediction advantage for the local
branch-adapted models. No labels were used to refit a model.

The mode plot adds +0.5% only for the two selected paths. Signed Qx/Qz keep
the original T Cartesian metric and ordered continuous lift. The midpoint
preserving shuffle remains near the zero-strain pattern; T keeps mixed
Qx/Qz with intermediate Qz. These observations describe the sampled
geometries rather than isolate one mode's causal contribution to energy.

The preserving zero-strain warm start already met the criterion at its
first target-plane evaluation; T took two optimizer updates. Finite
stopping/image resolution and path history remain limits on attributing
the observed nonlinear sampled response to continuous-path curvature.
The original .10 eV/A threshold and nine-image target remain as registered.

Replay from the repository root:

~~~sh
python -m scripts.analyze_hfo2_physical_holdout --output paper/VARNEB_JCTC/figures/hfo2_independent_20261018
~~~

This reads the frozen panel and completed exports; it does not contact HF,
run DFT or launch optimization. The public CSVs provide all90image records,
two outcomes, three selected-pair contrasts and all138forecast records.
All earlier training and preparation costs remain additional to the39
endpoint/path SCFs and48.808889 allocation core-hours of this physical test.

Figure5 compares complete paths, strain response and all numeric forecasts.
Figure3 overlays signed shuffles and open-axis stretches. Existing Python
matplotlib style is retained: normal-weight panel labels, four-sided frames,
aligned axes, large text and framed translucent white legends outside data.
Strain colors apply in5(a,b); channel colors in5(c,d), as stated in the caption.
The original training-only exports remain in the adjacent dated folder.

Verification:91scoped scientific/plotting/related manuscript checks passed.
The LuaLaTeX reading PDF compiled to48pages with no overfull/underfull boxes
or missing characters. Both standalone figures and PDF pages14,16,17 were
visually reviewed: captions stay with figures, axes/panel labels align,
legends do not cover data and the signed-error table remains readable.
This is not a full repository regression run or a new DFT calculation.
