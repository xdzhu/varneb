# Figure contract: provisional mechanisms, not final competing barriers

Core conclusion: same-PO+ snapshots distinguish a Pbcn-centred preserving
candidate from a Pbca-centred reversing candidate, while their scaled-cell
residuals and geometric shuffles show why ordinary convergence and full
mechanism/stability tests are still needed. Vanishing T-triplet amplitude
does not make a structure cubic. These are finite, unconverged observations.

Archetype: quantitative2x2 grid, with energy evidence(a) and verification/
mechanism panels(b–d). Backend: the existing Python/matplotlib workflow only.
This is a dated JCTC working-draft figure, not a submission-complete claim.

- (a) The three low-energy profiles in the commonPO+ thermodynamic view;
  ordinary-convergedPO→T, PO→Mstep20, preservingflipstep14. Show the actual
  relaxedT energy separately as a dashed reference, not a guessed saddle.
  The much higher reversingstep10 profile is tabulated and retained in full
  source data, not silently truncated or scaled to fit the low-energy panel.
- (b) All four source ordinary-NEB maximum-vector residuals per moving image
  with the unchanged0.10eV/Angstrom target. Fixed endpoints have no NEB force.
- (c) All four geometric rotatedT-patternQx traces. This is the shuffle sign
  relevant to the registered preserving/reversing variants, **not** polarization,
  a phonon occupation or an energy contribution. Allthree patterns are archived.
- (d) Frobenius norm of Green strain relative to the original orderedT cell.
  This is a dimensionless geometry measure, not an elastic energy or stress.

Final canvas183mm wide, approximately168mm high; PDF/SVG editable text,
PNG preview and compressed600dpiTIFF. User style takes precedence over
skill defaults:10pt ticks/legends,11pt axis labels, normal-weight13pt(a–d),
all four spines/ticks, no panel titles, aligned axes/labels, a shared boxed
semi-transparent white legend outside data. Straight connections only,
original discrete markers, no invented interpolation or error bars.

Source data:37evaluated image records, including reused fixed endpoints,
not37newDFTs or independent replicates. Freeze rawE/F/stress/INPUT/geometry/
logs and verify SinglePoint force replay; retain source direction/step/job.
Statistics: one deterministic snapshot per candidate, no mean/SD/CI or test.
Normalised arc uses each source's original joint metric; reversing a view
does not relift/reorder atoms or recompute the sourceNEB metric.

Review risks: only one source passes ordinary residual; no final barrier
ranking/strain-window/TS/phase-minimum claim. Three symmetry tolerances and
software warnings are reported, not tuned to preferred phase labels.
The Pbcncenter is currently a sampled energy maximum, unlike the stable
intermediate reported in the separately auditedSCAN study. Same space-group
symbol is not evidence of identical variants, pathway, stability or energetics.

QA: replay source fields before plotting; check shared legend/titles/spines,
panel/label alignments and text bounding boxes; inspect actualPNG in Python
output and confirm SVG retains text. Hash every exported artifact and CSV.
