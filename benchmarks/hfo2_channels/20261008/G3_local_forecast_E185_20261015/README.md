# E185: prospective local-face controls, not independent validation

2026-10-15. This freezes eight **uncomputed target Q values per original
center**, not a new material, a third raw Hessian center, a new strain, or
the physical +0.5% path holdout. Q denotes the existing two local coupled
coordinates (length in A), not a phonon wave vector.

## What is locked

Both T and preserving use their audited central nine training targets
Q1,Q2 in {-0.1,0,0.1} A. All nine labels exist at both centers; the outer
preserving holes are not filled. These labels and the earlier exploratory
comparisons were seen before selecting this test. Future labels were not.

The eight nonzero points of {-0.05,0,0.05} A squared are reserved as the
original eight independent local-face slots per center. No test point has
been computed. All five standard arms are retained: bilinear interpolation,
tensor quadratic interpolation, a refitted six-parameter total quadratic,
and both original full H/g quadratics (0.01 and 0.02 A probe steps).
The latter are the complete-space reference control, not a claimed new
harmonic theory. This is **80 scalar predictions, not 80 DFT evaluations**.

Each point has one precomputed safe start. Nearest central training Q is
selected by integer squared distance, then original ID, with the original
geometry guard. All 16 points select the known Q0 candidate. No new test
result becomes another test's parent. This tests that single-start branch;
it does not establish the global conditional minimum or branch uniqueness.

The manifest records raw energy changes relative to each raw peak,
origin biases, predictions calibrated using the already-visible released
Q0 energy, and origin-relative changes. Scoring retains signed errors,
maximum absolute error, RMSE, both steps and coverage out of eight.
The pre-label reporting target is **1 meV/f.u. = 4 meV/12-atom cell**.
This is an explicit engineering reporting target chosen after training,
not a numerical error estimate, physical guarantee or literature standard.
Unresolved differences cannot establish a predictive advantage.

## Fixed computation and admission

The original ABACUS/PBE/100 Ry/full 10 au DZP/Gamma 2x2x2/P=0/E=0,
zero-strain mechanical boundary and chart are inherited by exact source
pins. Original .5 A total displacement and 1.6 A distance guards remain.
Safeguarded BFGS/.03 eV/A released-gradient target/.05 A trust step/200
iterations are fixed; ordinary NEB .10 is unchanged. At most 12 **new
attempts**, including failures, per point: 192 for the two faces, no retry
or rollover. T34/81 and preserving35/81 modeling commitments are unchanged;
test slots reserved8/face, computed0/face. No allocation was submitted.

Before test labels: complete the existing B3/native released-Q0 curvature
audits and lock final candidate-model predictions. Both raw released R34
steps must be positive as necessary local admission evidence; this is not
a global stability or total-error certificate. Failure preserves the
uncomputed test and original prediction record; it does not authorize
changing its coordinates, electronic settings or branch to make it pass.
Two existing committed jobs still occupy the study allocation allowance.

This is **not** the full B0--B5/strong T-Cmma/physical-strain prediction or
final representative/domain freeze. These elementary energy controls
supplement, not replace, those comparisons. The +0.5% geometry and labels
remain sealed. Local-face evidence cannot substitute for path-barrier tests.

## Current numerical implication and cost

On the uncomputed points, maximum tensor-versus-H(0.01) origin-relative
prediction differences are T **0.502421**, preserving **1.465334 meV/cell**.
Those are model disagreements, not prediction errors. The test may not
resolve a winner; coordinates are retained rather than moved to obtain
an apparent gain. No independent improvement has been demonstrated.

The nine training targets historically cost 28/36 face SCFs (T/preserving),
not nine SCFs each. Full original faces cost 118/172, raw curvature296,
and Gamma/canary/G2/all branch preparation remain additional. Models have
different information and preparation costs. No equal-budget, DFT saving,
acceleration or B5 novelty claim is made.

## Reproduction and preserved amendment

`freeze.py --output <fresh-file.json>` reproduces `forecast_v2.json` byte
for byte. Original `forecast.json` (283f2b9d...) is retained unchanged.
The first independent checkout found 35 passes/1 metadata failure: existing
Windows CRLF source bytes differ from Git archive LF bytes. Revision2 only
uses canonical LF **text-code** hashes. All original numerical predictions,
starts and scoring/execution policies are checked unchanged. Native input,
DFT result, JSON and binary pins remain strict byte hashes. The raw failure
XML is retained; this was not a failed DFT or a change to scientific inputs.

Final related checks: 38 local and 38 independent, three modules, eleven new
checks; not the entire repository suite. The independent CLI produces
SHA256 aa77cdd7460ad809a5fbe06e80c684f7f73ce118835067c380ec219b96abd2d8.
Tested source tree04f584c738e66f3deebdaa42d61c089a7fe80248;
archive9ab42037c6320c7ad27ec1567c3f2902d10111fb59a1d5c0ce5a00230eeddff0.
Both snapshots and scratch are on C:, not the almost-full E: drive.

`status.json`: 08:55 CST both 28867177 and 28868827 remain quota-pending,
zero new labels. No new DFT/scheduler write, running-source overwrite,
private author-array access, parameter change, CI, release or installation.
User edits and manuscript are unchanged. Full research goal remains active
and incomplete; these forecasts are preparation for a falsifiable test.
