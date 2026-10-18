# Atom-resolved comparison, not a new mode assignment

`analysis.json` and `center_motion.csv` compare the original and actual E273
step9 ordinary band using the existing original-T proper rotation and
unchanged atom identities/periodic lift. No DFT, relabeling, mode fit or
synthetic atom interpolation was used to generate the nine samples.

Individual oxygen y displacements at the original center reach0.322A,
even though their Hf-relative species mean is nearly zero. The new/old
center difference includes oxygen z pairs near+/-0.200A and Hf x motion
near+/-0.068A. The mean transverse feature therefore does not describe
the whole rearrangement or establish a single eigenmode.

Each band's motions are measured from its unchanged initial endpoint in
the fixed T metric; these are not the absolute T-parent structural proxies
of Figure6. The small nonzero initial proxy explains that reference offset.
The table's lower-minus-original differences remove common translation;
Hf-relative columns instead remove mean Hf displacement.

Positive ordered fractional-z margins place each oxygen between the same
Hf layer bands at all nine samples and throughout their piecewise-linear
fractional segments. These are not physical distances or proof about an
unknown smooth MEP. No complete C:N/S:T family, phonon, electronic
polarization, TS index or energy partition is assigned. The descriptive
old/new comparison does not replace the matched E272 F/A/J intervention.

```text
python -m scripts.analyze_hfo2_oxygen_motion --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_oxygen_motion.py
```
