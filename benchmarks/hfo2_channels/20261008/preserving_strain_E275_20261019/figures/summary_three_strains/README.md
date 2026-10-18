# Internal accommodation changes the sampled strain response

This is manuscript Figure9, derived from the complete original two-target
outcome and the retained zero-strain/registered reference bands. It answers
one material question: does admitting the observed mixed rearrangement
change preserving reversal relative to escape under the same clamp?

The answer at the declared ordinary sampling level is yes. Released-source
peaks are83.166281/87.718545/91.698937meV/HfO2, rather than registered
96.828183/98.530877/125.424431. The0-to-+1% rise is8.532656 instead of
28.596248. T gives84.416624/90.418178/96.763650; the new small positive
T-minus-preserving gaps are not resolved physical barrier orderings.

(a) Native sampled-peak bars; (b) signed T-minus-preserving bars;
(c) transverse species-mean displacement at the actual peak;
(d) its open-y cell response. The nine rows in`source_data.csv` are actual
measurements, not fits, statistical replicates or extra strain points.
Categorical marker offsets in(c,d) are only for readability. No lines
connect or interpolate the strain response. Ordinary force tolerances
are not converted to energy error bars.

All nine inputs pass their own ordinary residual criterion. The new
midpoint needed zero updates; this is retained explicitly rather than
described as proof of an independently relaxed unique branch. A transverse
structural proxy is not Berry polarization or a named phonon. The failed
+1% local A/J intervention is not replaced by these free-chain results.

`qa.json` pins the complete analysis SHA256 and export contract. Python/
matplotlib produced editable PDF/SVG and300dpi PNG,183x155mm, normal13pt
panel labels, aligned axes/labels, top/right ticks, and an external framed
semi-transparent white legend. Rendered PNG was visually checked for
legibility, clipping and legend/data overlap.

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/plot_summary.py --source benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/completed_analysis_full/analysis.json --output FRESH_DIRECTORY
```
