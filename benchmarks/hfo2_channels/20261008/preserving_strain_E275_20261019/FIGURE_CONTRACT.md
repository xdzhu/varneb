# Same-strain mechanism comparison (existing Python workflow)

Question: does the released-center zero-strain candidate change the sampled
competition after transfer to an existing substrate boundary? The plot must
show the actual energy profile together with atomic and cell accommodation,
not merely the smallest scalar peak. This is a quantitative2x2grid.

- (a) Native energies from the same exact initial cache, each candidate's
  own normalized generalized arc. No spline, synthetic images or peak fit.
- (b) Fixed-original-T mean O-minus-Hf y displacement, a structural proxy,
  not electronic polarization or an assigned phonon eigenmode.
- (c) Open-axis response100(Fyy-1), using that same fixed reference.
- (d) Actual ordinary NEB per-interior residuals; endpoint residuals are
  undefined, not zero. Dotted line denotes the0.10target.

Include registered preserving, released-source preserving and T escape;
the independent upper-initialized preserving control is also shown at+.5%.
Never overwrite the old labels or pick only the lowest initializations.
Only an actually collected terminal strain can be drawn by this entry point.
The current+.5% target used7fresh SCFs, zero optimizer updates; its ordinary
force pass does not prove initializer-independent relaxation or a smooth MEP.
At+.5%,87.719vs90.418meV/fu is a close sampled competition, not a certified
2.700meV barrier ordering. +1% remains pending until its original handle ends.

Exports:183x155mm, editable PDF/SVG plus300dpi PNG, exact source CSV,
analysis SHA and Python QA JSON. Labels(a)--(d)normal13pt, larger tick/axis
fonts, top/right axes and inward ticks, aligned grid, no subplot titles.
Shared legend outside curves, semi-transparent white background and border.
One deterministic band per source; no replicates, error bars, statistics or
energy uncertainty inferred from the force threshold. No raster editing.
