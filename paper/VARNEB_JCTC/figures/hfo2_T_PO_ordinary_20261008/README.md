# Audited ordinary HfO2 T--PO band

Reproduce from the repository root with the existing Python/ASE environment:

```powershell
$env:MPLBACKEND = 'Agg'
$env:OMP_NUM_THREADS = '1'
python -m scripts.plot_hfo2_converged_band --output <new-directory>
```

The output must be a new directory. All ten frozen SCFs and every archived
descriptor are re-audited/replayed before plotting; no DFT is called. Straight
connections do not interpolate a saddle. The source CSV retains all six Green
strain components, while the figure shows their three diagonal components.
Blank CSV fields denote undefined zero-reference fractions or fixed-endpoint
NEB residuals, not measured zero values.

The script initially labels visual review pending. This checked bundle's QA
records the completed PNG inspection and separate vector-text checks:61 SVG
text elements and446 extractable PDF characters, with three embedded TrueType
fonts. All(a)--(f) labels are normal weight, axes/labels aligned, all spines
visible, no panel titles, and framed translucent white legends.

The figure contract and scientific caption are in the research manuscript
directory. This is one ordinary-residual-passed edge; full-variable TS,
sampling/error bounds and competing-channel predictions remain incomplete.
