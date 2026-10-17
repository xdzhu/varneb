# E265: one finite, zero-new-DFT anharmonic model experiment

Registered after all E260 training labels and before every remaining P8 label.
This is model selection on seen training data, not a pre-training prediction.
It adds no training candidates, DFT calls, Hessian centers or independent
points. Original protocols, forecasts, refusals and all preparation costs stay.

Use all81 E260 native energy/full-gradient records, including rejected trials
and capped candidates, as **physical point observations**, never relabelled
as conditional minima. Pin the E264 complete native report. The exact origin
is the old released preserving candidate; its energy/full gradient are native.
Coordinates are displacement from that origin in the original36D length chart.
Retain the two original Q axes and the exact frozen zeta direction. The other
33 directions are their orthogonal complement; they are not new phonon irreps.

For each original E254 step0.01/0.02A independently, project the full measured
34-column physical-gradient action into this chart. Apply the already declared
reciprocity correction/extension, preserving the raw data and defect. This
determines the free34 block and its cross terms with Q. The unmeasured symmetric
Q/Q block has three **fitted** coefficients, not an assumed measured zero.
Do not combine curvature columns from different centers or clip eigenvalues.
The raw R33 restriction must be positive; its value is a local source-point
property, not certification throughout the fitted domain.

Compare exactly two families per step, no degree/hyperparameter search:

1. Harmonic full-space reference: fit only the three Q/Q coefficients.
2. Same base plus a retained-direction nonlinear residual: zeta^3,zeta^4,
   Q_i*zeta^2,Q_i*Q_j*zeta^2 (seven terms), and for each released coordinate r_j,
   r_j*zeta^2,r_j*Q1*zeta^2,r_j*Q2*zeta^2 (99 terms). Total109 coefficients,
   versus3 for the harmonic reference. E_origin and g_origin are fixed inputs.

This is an ordinary gradient-enhanced polynomial/harmonic-bath ansatz, not a
new theorem or a demonstrated advantage. It deliberately tests whether the
observed coupled relaxation is representable with a bounded nonlinearity.
Complete fixed-reference re-expression cannot produce an independent gain;
the full T-optical+cell harmonic reference remains the same full-space solve.
All original restricted T/Cmma controls and both steps must remain reported.

Normalize coordinates by0.1A. Least squares uses energy and all36 physical
gradient components, multiplying gradient rows by0.1A to match energy units.
Normalize design columns and use SVD relative cutoff1e-10. Rank-deficient
fits are refused, not silently used as minimum-norm forecasts. Report rank,
conditioning, energy and full-gradient residuals, not only energy agreement.
Report five leave-one-Q-family-out **training** validations, including the
central/capped family; do not count them as new independent data or choose
the winning step after seeing P8. Neither training error nor1meV/fu reporting
target is a physical DFT-error bound.

The experiment first answers identifiability and training support. It does
not automatically freeze a final model, extrapolate across the failed-center
hole as a certified stable sheet, certify branch minima/TS, or authorize P8
or physical+0.005 DFT. Model/domain/strong-control/test-target and prospective
predictions must still be locked before any remaining independent label.
If this one experiment fails, retain that finite result; no automatic new
degree, radius, extra training, full Hessian or material matrix follows.
