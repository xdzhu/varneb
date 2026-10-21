# HfO2 final two-cross closure

User approved the finite blueprint on 2026-10-21. Active goal: finish the
two missing strain-pattern crosses, interpret either outcome, integrate the
manuscript/main figure/case, compile and show the reading PDF, then push.
No further material, strain, model or Hessian matrix is authorized.

## Fixed scientific argument and terminology

Question: how does internal atomic accommodation under partial clamping
change homogeneous reversal relative to escape from the same polar well?
Existing whole-path evidence reduces the sampled 0-to-1% preserving penalty
from 28.596 to 8.533 meV/HfO2. Local fixed-pattern controls explain part,
but not all, of that change. The final experiment separates imposed strain
from changing the four retained geometric observables.

- Preserving reversal: preserves the signed reference shuffle, not every
  oxygen atom's coordination identity.
- q0/q1: four original-T geometric observables at the original zero/+1%
  sampled peaks; not the earlier two local joint coordinates Q1/Q2.
- J: joint release of the 29 remaining atomic and three permitted cell
  directions, with three translations gauge-fixed.
- Transverse proxy: mean O-minus-Hf displacement, not electronic polarization
  or an independently identified phonon eigenmode.
- U: local candidate energy above the same-strain initial well, not a
  newly calculated path barrier or a certified saddle.

## Finite execution

Protocol/case: `benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021`.
Two existing minus-sign J diagonal candidates are reused. New conditions:
zero strain with q1 and +1% with q0. Hold all four values exactly.
Electronic inputs unchanged: original 100 Ry/full10auDZP/PBE/six assets/
Gamma2x2x2/zero pressure and field. At most 40 SCFs per cross, 80 total;
two one-node32MPI hfacnormal01 allocations, 2h each. No new A treatment,
second sign, midpoint, whole path or automatic numerical-budget extension.

Source `1a4925da` was pushed before allocation release. Bundle pin:
`ecb498b40ec3f19b593bfca198a754d1f2cdc62deb389fcc1dbb6b0d85ea2683`.
Remote namespace:
`/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1`.
Zero-SCF native preflight passed both diagonal caches and cross geometries.
Original held jobs: 28995326 (zero_q1), 28995327 (plus1_q0).
Both held-resource checks passed and the original jobs were released.
Slurm reports both RUNNING from 12:52:27 CST, on node4 and node6 respectively.
This is allocation startup, not yet evidence of completed electronic calls.
Follow only these two original handles; do not resubmit.

## Paper closure and stop rule

Use research-paper, English/JCTC framing. Evidence-first update order:
Results -> Discussion/Conclusions -> Introduction/Abstract. Keep the four
approved Results sections and the existing evidence tables/figures, with
old negative prediction results identified as completed rather than an
unfinished promise. Do not ask the author to reconfirm the approved scope.

The final cross figure uses the existing Python/matplotlib workflow, source
data only, normal-weight (a) panel labels, top/right axes and large fonts.
Four panels link local energy, pattern contrast, transverse displacement
and open-axis response. No mock data or fitted surface is required.

Accept terminal results whether positive, negative or inconclusive; no
automatic expansion. Then update the concise Chinese reading guide and
reproducible case, run relevant native/constraint/figure/manuscript tests,
compile with the existing LaTeX installation, visually inspect changed
pages and present the PDF. Push the completed evidence and source.
Preserve the five unrelated pre-existing user edits in the working tree.
