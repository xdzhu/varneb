# E184 — elementary nonlinear controls on already seen candidate sheets

Previous E183 was PROGRESS/28bc8b29. Question: does an unchanged peak-local
quadratic remain a competitive finite-domain energy model once elementary
energy-only comparators see coarse candidate energies? This is an
**exploratory post-hoc coarsening experiment**, not prospective validation,
new DFT, a new VCNEB algorithm or proof of H2.

## Fixed split and comparisons

Use the original Q grid(-.2,-.1,0,.1,.2)Angstrom. Both integer indices even
select nine coarse coordinates; the remaining sixteen are fine coordinates.
Selection uses geometry only, not error. All labels were already seen before
this experiment; these Q use the original joint length chart, not phonon
mass amplitudes. Model families were chosen after prior training results
were seen; coefficients use coarse labels only. All tested arms are retained:
bilinear interpolation, tensor-product degree2in each
axis interpolation, a six-coefficient total-degree2least-squares fit, and
the unchanged H/g models at .01/.02probe amplitudes. There are five reported
model arms because the original quadratic has two probe steps.

Tensor degree2includes mixed terms up to total degree4; it is not another
name for a harmonic quadratic. These are standard comparators, not new
theory. Interpolation requires every nonzero-weight coarse label. Exact
zero weights need no unavailable label; negative Lagrange weights remain
visible. Rank-deficient global quadratic fits abstain, with no unidentified
coefficient/pseudoinverse fallback. No forces are treated as exact envelope
derivatives, no geometry is repaired and no missing point is filled.

## Actual T result on the same sixteen fine coordinates

Origin-aligned errors below use **meV per twelve-atom Hf4O8cell**, not
meV/HfO2 and not a physical/whole-path barrier error. Each model subtracts
its own seen-origin prediction and the observations subtract their seen
origin; the report also retains signed unaligned errors and origin biases.

| Model | Max absolute shape discrepancy | Shape RMSE |
|---|---:|---:|
| Bilinear | 6.900894 | 4.219064 |
| Tensor quadratic | 2.003468 | 1.043743 |
| Refitted total quadratic | 7.667090 | 4.006016 |
| Original H, .01 | 12.901228 | 5.177484 |
| Original H, .02 | 12.975249 | 5.219402 |

The refitted quadratic has a +1.546418meV/cell origin bias; its unaligned
maxerror6.120672differs from its shape maximum. Both are retained. No model
or step is discarded because its result is less favorable. Lower tensor
error shows that elementary finite-coordinate interpolation is a relevant
control here, not that a unique anharmonic mechanism has been identified.
Interpolation versus the unfitted peak-local H also changes information
and preparation, so this is not an equal-information advantage of B5.

## Preserving coverage is not a two-dimensional success

Only five of nine coarse labels and twelve of sixteen fine labels exist.
The global six-coefficient quadratic design has rank5and is not fitted.
Consequently the **all-five-model common set is empty**; statistics are null,
not zero. Without that unidentifiable fit, the four other arms share only
four axial coordinates, ids1–4. The bilinear/tensor maxima there are
0.812719/0.470404meV/cell, versus4.671404/4.720758for the original two H models.
On these axes, tensor interpolation reduces to a one-dimensional quadratic;
it does not validate a two-dimensional coupled landscape. The other eight
observed fine labels and four fine holes are retained, not quietly removed
from the coverage record. No claim that every quadratic model fails follows
from the missing global fit.

## Why this is not a cheap independent prediction

All eight non-origin coarse starts in **both** channels inherit a fine-node
continuation parent. Holding fine energies out of the fit therefore does not
blind the coarse geometries, branch history or preparation to fine data.
This is a functional interpolation diagnostic, not an independent acquisition
or clean cross-validation of the whole physical workflow.

Historical calls at the selected coarse/fine coordinates are47/71forT and
58/114forpreserving, totaling the original118/172SCFs. Preserving coarse
costs include46calls at four missing-label nodes; fine holes consumed48,
94failed-to-converge-node calls total. A missing coarse corner stopped after
ten, not the twelve-attempt ceiling. Original296curvature/Gamma/canary/G2/
seed preparation remain additional. Nine coarse coordinates are not nine
SCFs and these historical subsets are not a prospective standalone cost.
There is no demonstrated DFT saving, whole-band acceleration or physical
accuracy bound. The original .03released-gradient criterion is not changed.

Decision: retain elementary local energy interpolation as an additional
comparator before attributing future independent gains to an adaptive model.
It supplements, not replaces, B0/B1strain-barrier or strong T/Cmma controls.
No predictive domain/model/threshold for the original independent tests is
frozen by these exploratory results. Stable branches, measured resolution
and independent face/path predictions remain required.

## Reproducibility and state

E180scalar source39e31d79and complete native face proofs are hash-pinned;
no raw-SCF/private author reparse or new reference selection.37local/
37independent related checks over3modules,12new,0finalfail/errors/skips,
notfullsuite. Independent CLI reproduces complete analysis bytes
377f2c3e2b6c21f0f5d033b2e180ed233fea7547924f2ac5ab993aa4b041f8a5.
Testedtree3bbd688b791371794eeaf5431a203a065685be91,
archiveca97e2dc90939bc215fc16d94506261a18628ab059dfd458dea9a0c05b562f83.
First36pass/1fail was an incorrect **test expectation** that coarse holes
used48calls; source-derived arithmetic already used the actual46. Corrected
the expectation from the pinned native proof, not the energies/analysis.
Raw first failure is retained; the initial report before detailed cost
fields remains in Cscratch. All new fixtures/archives use explicit Cscratch.

08:36:18CST28867177/28868827stillPENDING/AssocGrpCpuLimit,0newSCFs/labels.
No new calculation, reservation, scheduling write, material, path, parameter,
radius, threshold, CI/install/release/PyPI/monitor/submission. Independent
eight face slots each and physical+.005remain sealed. Main92f15239and five
user edits preserved. Full goal active/incomplete; no new certified PES/TS.

```text
python benchmarks/hfo2_channels/20261008/G3_coarsening_E184_20261015/analyze.py --output fresh_report.json
python -m pytest -q tests/test_hfo2_coarsening.py
```
