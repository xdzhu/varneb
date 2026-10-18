# Same +0.5% substrate, different sampled preserving candidates

The released-source candidate has87.718545meV/HfO2 peak versus98.530877
for the registered preserving band and109.351721 for the upper-initialized
preserving control. The T escape remains90.418178. All four use the same
exact native initial well; all three preserving candidates also use the
same exact final well. These are deterministic sampled profiles, not replicates.

Panels(a)--(d) show native energies, the original-T transverse structural
proxy, open-axis response100(Fyy-1) and per-interior ordinary residuals.
Each line uses its own native generalized arc. Dotted line in(d) marks.10;
endpoint residuals are undefined and excluded. Straight connectors are not
a smooth-MEP interpolation. The structural proxy is not electronic polarization.

The new candidate used7fresh interior SCFs and zero optimizer updates;
residual.092315 passed the same ordinary criterion. Its2.699634meV/fu
separation from T is not a certified physical ordering or energy error bound.
This is consistent with the original candidate comparison changing when
mixed accommodation is admitted, not proof of a unique branch or TS.
+1% is still absent from this partial-terminal comparison, not assigned zero.

Python/matplotlib generated editable PDF/SVG and300dpi PNG, all from
`completed_analysis_0050/analysis.json`; its SHA is recorded by`qa.json`.
The36rows in`source_data.csv` are the exact plotted samples. Normal13pt
panel labels, top/right axes, aligned labels and an external framed white
semi-transparent legend follow the author contract. Visual review found
no clipping, legend/data obstruction or unreadable labels.

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/plot_comparison.py --source benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/completed_analysis_0050/analysis.json --strain .005 --output FRESH_DIRECTORY
```
