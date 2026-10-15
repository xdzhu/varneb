# E188: all registered reference controls at the locked independent targets

2026-10-15. E187 was a verified wait, not new science or a code milestone.
This turn freezes previously missing **reference-control predictions** for
the E185 test coordinates. It does not generate/read a test DFT label.

## Fixed information and scope

Two original centers x two original derivative amplitudes x seven references
x eight already locked targets =224 scalar model-point records. Reference
inventory is unchanged: T geometric triplet(rank3), two complete optical
doublets(rank4), all33 T optical directions, and Cmma frames103/106/114/127
(each rank4). Each includes the same three open-cell directions. All use
the same measured PBE H/g, shared Q constraints, mechanical boundary and
visible released-Q0 energy calibration. No fitted coefficient, negative-mode
clipping, pseudoinverse, deleted Cmma direction or post-label frame selection.

These are the E180 **restricted shared-Q** controls: omitted directions are
clamped. They are not a full outside-reference release or a completed native
Cmma-coordinate conditional model. Author LDA mode vectors are geometric
directions only, not transplanted PBE frequencies/curvatures. Private raw and
mapped mode/response arrays remain in memory and are not redistributed.

All raw formal energy changes, origin biases and calibrated predictions are
retained. A target outside the original geometry guard has **null accepted
prediction**, not an extrapolation scored as a poor DFT prediction. The
original E185 target coordinates, starts, 1meV/f.u. engineering reporting
target,192-attempt ceiling and all scoring/coverage rules remain unchanged.

## What the actual calculation found

At both probe amplitudes, geometry-eligible target counts out of eight are:

| Center | T triplet | T two doublets | T all optical | Each of four Cmma frames |
|---|---:|---:|---:|---:|
| T | 6 | 8 | 8 | 6 |
| preserving | 0 | 8 | 8 | 2 |

Thus140/224 accepted scalar predictions,84 geometry-domain abstentions;
**zero new DFT points**. All calibration origins pass the geometry guard.
The all-reference common set on T is only the four axial targets1--4, not
a 2D validation. On preserving it is empty. Each Cmma frame there supports
only x-axis targets1/4; the triplet requires formal displacements of45--98A
to satisfy some shared Q constraints and is rejected. This is a limitation
of that restricted span for this Q/domain, not proof that every Cmma/T native
model fails. The full T and two-doublet controls both cover all eight targets
including off-axis points and must not be omitted from future comparisons.

Full T33+cell3 and the already-frozen full-space Schur/H model agree within
1.6653e-16eV/cell: standard same-problem algebra, **not predictive gain**.
They are not independent algorithm replications. Likewise, four Cmma gauge
registrations are robustness alternatives, not statistical repeat samples.

On the same eight accepted targets, two-doublet vs full-T calibrated model
disagreements are474.975/474.958meV/cell(T) and542.960/542.668meV/cell
(preserving), at .01/.02A steps respectively. Those are mathematical model
disagreements, not measured DFT errors or physical barriers. Different
allowed spans and clamped directions change the constrained response.
Accuracy/advantage still needs the independent labels and numerical resolution.
All energies are per twelve-atom Hf4O8 cell; divide by4 for eV/HfO2.

## Verification, cost and next gate

49 local and49 independent related checks, four modules/nine new, zero final
failures/errors/skips; not the full suite. All224 formal energies and squared
coordinate norms were separately reconstructed from the original E180 scalar
polynomials, without private vectors in the tests. The independent CLI with
read-only private sources reproduces forecast.json bytewise:
97692cb0989a4bafdff86349a732126cce5d7d93d9a7d47926ce844fb1f0b285.
Tested source treefd0871fafeff181d5b8a727aac5cd4eb6dbe80a5;
archive84b6704b9a81713d1f85920c502b83574eda2c574e209cf9732b9db04856b347.
All scratch/snapshots are on C:; no writes to the private E: source bank.

Reproduce with `freeze.py --source-root <private-author-bank> --output <fresh-file>`.
Existing outputs and output paths inside the private bank are refused before
reading private sources. Source/prediction/code pins are in forecast.json.

Original296 curvature SCFs,118/172 full training-face calls, Gamma, canary,
G2, reference mapping and branch preparation remain additional. Same measured
H does not imply equal preparation cost; no DFT saving/acceleration/novel
harmonic theory claim. Eight independent slots/face remain reserved/computed0,
training T34/81 and preserving35/81 unchanged, physical +0.5% still sealed.

09:20:20 CST:28867177/28868827 are still live quota-pending,0newSCFs. Existing
B3/released-Q0 raw curvature audits and final candidate-model predictions
are still required; this is not the full B0--B5/physical-strain/native-chart
freeze, production admission, stable conditional branch or JCTC acceptance.
No new DFT, allocation, scheduler write, physical parameter, live runtime,
CI/install/release/monitor or manuscript change. Five user edits preserved.
The bounded experiment workflow keeps coverage and abstentions explicit,
instead of selecting a favorable reference or silently replacing forecasts.
