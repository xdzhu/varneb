# E080 paired clamped training-band figure contract

Core conclusion: geometrical displacement coverage and atomic energy-gradient
coverage answer different questions, and neither is transferable between escape
channels without testing. No energy partition or predictive superiority claim.

Evidence chain: actual terminal T audit E079 and M audit E068 -> nine native
portable E/F/stress records per band -> fixed E023 T/Cmma references and
inverse prescribed G2 frame -> E080 T and historical E071/E078 M scalar reports
-> pinned CSV -> Python figure -> visual/structural QA -> main manuscript.
Native SCF/ordinary-residual replay happens in analysis, not in the plotter.

Archetype: 3 rows x 2 columns. Left PO+ to T, right PO+ to M. Top: nine
discrete P=0 energy samples; middle: squared common-T displacement fractions;
bottom: squared translation-free atomic-gradient fractions. All four Cmma
registrations remain. Bottom and middle share 0--100% scale but have different
denominators. Each column uses its own normalized original generalized arc;
no claim that equal normalized positions are the same physical structure.

Backend: existing repository Python/Matplotlib scientific plotting workflow;
no R switch, AI image generation, DFT launch, fit, smoothing or new mode choice.
Export vector PDF/SVG, review PNG300dpi, lossless TIFF600dpi and scalar CSV/QA.
Paper uses PDF. Keep historical M-only E072 assets unchanged.

User style: normal-weight large (a)--(f); top/right spines and ticks;
aligned equal-size axes in each grid row/column; no subplot titles;
legends outside curves with translucent white background and visible frame;
readable >=10pt ticks/labels and >=11pt labels. Shared external reference
legend, separate channel labels in top-panel legends, sampled peak guide only.

Risks and controls:

- T endpoint has exactly zero displacement from T: its fraction is undefined,
  not zero physical coverage. Keep the source zero-denominator flag in CSV;
  omit that endpoint from displacement curves instead of fabricating a value.
- Four Cmma-frame values are registrations, not independent replicates or
  uncertainty estimates. Draw all four; no statistical confidence band.
- Rank3 and rank4 are not equal-dimensional accuracy controls. Labels retain
  ranks. No unmeasured cell gradient is folded into atomic fractions.
- Near-zero gradients can make percentages sensitive. CSV retains absolute
  internal/retained/omitted norms, zero-gradient and translation diagnostics;
  fractions are descriptive, never a force or energy error certificate.
- Source phonons/LDA directions are not new PBE eigenpairs or local TS modes.
- Same original inputs, mechanical boundary and cached ordered PO+ record
  must be checked; reference/receipt/content pins reject changed source data.
- Ordinary residual pass and a highest sampled image do not certify a saddle,
  conditional branch or a barrier ordering. No G3 selection or holdout reads.

QA acceptance: real scalar values/units and pins verified; zero endpoint
masked; all registrations retained; shared limits/alignment/label styles
tested; fresh-only exports; actual PNG and manuscript placement inspected.
