# Complete matched atomic-release experiment

2026-10-15. The original job28867177 finished normally on node6/hfacnormal01:
4577s, one node/32CPUs, all10atomic-gradient candidates,36new SCFs. The
unchanged E176 complete auditor reparsed all36native records, checked old
runtime/input provenance and exactly replayed every original optimizer
trace. This supersedes the *coverage*, not the preserved bytes or values,
of E192's seven-node prefix.

## What the experiment now tells us

对翻转路径附近原先指定的五个局部点，主要降能均已发生在原子重新排布
这一步：原子释放降能44.06--48.35meV/12原子胞，随后允许晶格与原子联合
放松，额外降能0.80--1.14meV/胞。因此，不能把这组局部点的大幅降能主要
归因于额外的晶格释放。局部二次模型在原子释放后也已出现明显的面形偏差，
不是直到放开晶格才出现。下一关仍是检查释放点的横向稳定性；这组结果
不是整条路径的机制结论、唯一模式能量分解或独立预测证据。

F freezes omitted coordinates, A releases31atomic directions at fixed Q
and the frozen-Q seed cell, and J jointly releases atomic and permitted
open-cell coordinates. Q itself can contain cell motion, so different Qs
do not define one common fixed cell. These are treatment differences between
local candidates, not uniquely separable atomic and strain energies.

All following energies are **meV per12atomHf4O8cell**, four formula units.
Divide by four for meV/HfO2. Ordinary NEB0.10eV/A is unchanged; the separate
released-gradient criterion is0.03eV/A. No electronic input was retuned.

| Preserving Q1,Q2/A | A-F | J-A |
|---|---:|---:|
| 0,0 | -48.348691 | -1.136239 |
| -0.1,0 | -44.073849 | -0.804122 |
| 0,-0.1 | -47.241146 | -1.006188 |
| 0,+0.1 | -47.238928 | -1.006405 |
| +0.1,0 | -44.063824 | -0.816911 |

At the T origin, A-F=-0.271398 and J-A=-0.380744meV/cell, a much smaller
local release. The original twenty model-point predictions, both finite
probe steps, origins and signed errors are retained without refitting.

## Matched four-axial-point model comparison

The following are maximum discrepancies of *origin-aligned energy changes*
on the same four nonzero axial points at each center, not independent test
errors or measured physical uncertainty:

| Center | Original Hessian step/A | F | A | J |
|---|---:|---:|---:|---:|
| T | 0.01 | 0.416264 | 0.501879 | 1.327177 |
| T | 0.02 | 0.403160 | 0.515704 | 1.346071 |
| Preserving | 0.01 | 0.075880 | 4.349824 | 4.671404 |
| Preserving | 0.02 | 0.116535 | 4.399770 | 4.720758 |

For preserving, the largest missed origin-aligned A-F treatment change is
4.284867/4.284892meV/cell for the two models; the corresponding J-A
discrepancy is0.332117/0.331406. The original postprocessor checks the
signed **pointwise** decomposition, not a subtraction of maxima from
different points. These findings localize the measured treatment mismatch
but do not uniquely distinguish anharmonic atomic response, finite residuals
and branch changes. No step is selected as the winner and no coefficient
is fitted after labels.

All10atomic candidates have joint released-gradient norms above0.03eV/A;
none becomes jointly stationary just by passing the atomic-only test.
The released-Q0 curvature job remains a separate unfinished prerequisite.
No global conditional minimum, stable surface, index-one TS, complete
B0--B5 freeze, independent gain or JCTC readiness is claimed.

## Evidence, costs and reproduction

- atomic_native_complete.json:8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92.
- analysis.json:3cc150ec3a6d0dbfa5f5879f7215a3692524309add3e5bb7f5abed384a197444.
- Original E176 complete reader and E178 analysis code are unchanged.

```text
python benchmarks/hfo2_channels/20261008/G3_treatment_comparison_E178_20261015/analyze.py --atomic-audit benchmarks/hfo2_channels/20261008/G3_atomic_complete_E194_20261015/atomic_native_complete.json --atomic-audit-sha256 8a2fdf2a2beec306bb2f39aed8b347019a4f049b2d8cbc42cd79de1fce922b92 --output NEW_COMPLETE_COMPARISON.json
```

New atomic attempts/SCFs36:9atT and27atpreserving. Measured completed-call
proxy39.749558core-hours; original allocation40.684444core-hours. Frozen
statics8/8.136814proxy core-hours, original curvature296SCFs, previous joint
paths/canaries/preparation and the new272-probe curvature batch are additional.
Neither budget120versus observed36 nor comparison with a censored failed
optimizer is a speedup measurement. No total campaign-cost claim is made.

The ABACUS workflow skill guided original-asset/native EFS/MPI provenance
validation. All execution stayed on the authorized HF workflow, not generic
235/direct-node examples. No new DFT or scheduler operation was initiated,
production source changed, parameter/guard/budget adjusted, independent
slot consumed, private reference array read or physical+.005holdout opened.
Main manuscript and five user edits remain unchanged. The full goal is active.

Related verification:83local and83independent checks across five modules,
including eight new actual-complete-cohort checks, zero failures/errors/skips;
not the entire repository suite. The independent CLI reproduces analysis.json
bytewise. Tested source treedc0e7f0538ad30b177b7ca20d69a6d049db82acf/archive
08b0ebfcdf2fd6dc9af0772c8a9be631a257aea4b78cc0de0f2554094583e09c.
All scratch and independent exports are on C:. At11:08CST,28868827 remains
RUNNINGnode9/32CPUs,41/272SCFs and no failure record. Its full derivative
results and the final independent-test gates remain missing. No new test
allocation was admitted just because the B3 prerequisite finished.
